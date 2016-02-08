#!/usr/bin/python

"""first check to see if the csv and all files are present
process the csv file to remove - . / and spaces
make sure the csv file has the correct number of columns
run the bcl script
wait
check to see if the nohup outputs the correct ending
email that done
"""

import os
import sys
import string
import random
import subprocess
import glob
import logging
from datetime import datetime
import smtplib
import ConfigParser
import time
import re

def main():
    # Open config file.
    config = ConfigParser.SafeConfigParser()
    if len(sys.argv) == 2:
        config.readfp(open(sys.argv[1]))
    else:
        config.readfp(open('pathway.cfg'))

    # Create logger.
    logger = logging.getLogger(sys.argv[0])
    fh = logging.FileHandler(config.get('Globals', 'LogFile'))
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    logger.info('Starting %s...' % sys.argv[0])
 
#find eligible runs logic starts here 

    input_seq_dir_name = config.get('Globals', 'InputSeqDirectory')
    already_run_f_name = config.get('Globals', 'AlreadyRunList')

    try:
        with open(already_run_f_name) as f:
            already_run_dirs = [z.strip() for z in f.readlines()]
    except IOError:
        logger.error(
            'Cannot find already-run file. Check AlreadyRunList in '
            'pathway.cfg.'
        )
        raise
    try:
        input_seq_dirs = os.listdir(input_seq_dir_name)
    except OSError:
        logger.error(
            'Cannot find read directory. Check InputSeqDirectory in '
            'pathway.cfg.'
        )
        raise

    # For each input directory, for each index length,
    # create sample sheet and run bcl2fastq.
    bcl2fastq_output_dirs = []
    for d in input_seq_dirs:
        d = os.path.join(input_seq_dir_name, d)
        if d in already_run_dirs:
            #logger.info('Skipping directory %s. (Has already been run.)' % d)
            continue
        logger.info('Processing %s...' % d)

        # Check preconditions.
        if config.get('find_eligible_runs', 'Locked') == 'True':
            logger.warn('Lock set. Not able to process %s.' % d)
            continue
        config.set('find_eligible_runs', 'Locked', 'True')

        looks_like_illumina_dir = has_correct_illumina_dir_form(d, logger)
        has_req_files = has_required_files(
            d, config.get('find_eligible_runs', 'CheckFiles'), logger
        )
        sample_sheet_path = check_for_single_csv(d, logger)

        if not looks_like_illumina_dir:
            logger.info(
                'Adding %s to already-run file (%s).' %
                (d, already_run_f_name)
            )
            with open(already_run_f_name, 'a') as f:
                f.write('%s\n' % d)
        if not (
            looks_like_illumina_dir and has_req_files and sample_sheet_path
        ):
            logger.info('Skipping %s.' % d)
            config.set('find_eligible_runs', 'Locked', 'False')
            continue

        # Run bcl2fastq for each index length.
        proc_csv_paths = []
        index_lengths_present = []
        for s in [
            int(z) for z in config.get(
                'find_eligible_runs', 'IndexSizes'
            ).split(',')
        ]:
            # Create sample sheet CSV file for index length.
            processed_csv_path = process_sample_sheet(
                s, d, sample_sheet_path, logger
            )
            if processed_csv_path:
                proc_csv_paths.append(processed_csv_path)
                index_lengths_present.append(s)

        dual_index_present = any([z > 8 for z in index_lengths_present])
        for i, s in enumerate(index_lengths_present):
            processed_csv_path = proc_csv_paths[i]
            # Run bcl2fastq.
            unaligned_dir = run_bcl2fastq(
                d, s, processed_csv_path, logger, dual_index_present
            )
            output_dirs_f_name = config.get(
                'find_eligible_runs', 'OutputDirsList'
            )
            logger.info(
                'Adding %s to list of bcl2fastq output directories (%s).' %
                (unaligned_dir, output_dirs_f_name)
            )
            with open(output_dirs_f_name, 'a') as output_dirs_f:
                output_dirs_f.write('%s\n' % unaligned_dir)

            bcl2fastq_output_dirs.append(unaligned_dir)

            logger.info(
                'Writing directory %s to already-run file (%s).' %
                (d, already_run_f_name)
            )
            with open(already_run_f_name, 'a') as f:
                f.write('%s\n' % d)

        config.set('find_eligible_runs', 'Locked', 'False')

    if bcl2fastq_output_dirs:
        logger.info(
            'bcl2fastq output directories to be processed: %s.' %
            ', '.join(bcl2fastq_output_dirs)
        )
    else:
        logger.info('No directories to be processed.')
     
                  
## 'start makes' logic starts here

# Read in values from .cfg file (pathway.cfg unless otherwise specified).
    try:
        with open(config.get('find_eligible_runs', 'OutputDirsList')) as f:
            unaligned_dirs = [z.strip() for z in f.readlines()]
    except IOError:
        logger.error(
            'Cannot find file containing list of bcl2fastq output '
            'directories. Check OutputDirsList in pathway.cfg.'
        )
        raise
    logger.info('%s possible runs to process.' % (len(unaligned_dirs)))
    not_done_dirs = []
    seen_new_proj_dirs = []
    for d in unaligned_dirs:
        logger.info('Processing %s...' % d)
        

        if (
            not os.path.exists(os.path.join(d, 'rsync_is_finished'))
        ):
            open(os.path.join(d, 'rsync_is_running'), 'w').close()
            rsync_files(d, config, logger)
            open(os.path.join(d, 'rsync_is_finished'), 'w').close()
        if (
            not os.path.exists(os.path.join(d, 'linking_is_finished')) and
            os.path.exists(os.path.join(d, 'rsync_is_finished'))
        ):
            open(os.path.join(d, 'linking_is_in_progress'), 'w').close()
            make_links(d, config, logger, seen_new_proj_dirs)
            emailed = email_results(d, config, logger)
            if emailed:
                if os.path.islink(d):
                    os.remove(d)
                open(os.path.join(d, 'linking_is_finished'), 'w').close()

        if os.path.exists(os.path.join(d, 'linking_is_finished')):
            logger.info('Finished processing %s.' % d)
        else:
            not_done_dirs.append(d)
            logger.info('Not done processing %s.' % d)

    with open(config.get('find_eligible_runs', 'OutputDirsList'), 'w') as f:
        for d in not_done_dirs:
            f.write('%s\n' % d)     
                                             
def ssh_project(
    new_project_directory, investigator, id_number, new_name, config, logger
):
    """Connect to pageGenHost via ssh, and run script located at
    pageGenPath.
    """
    if os.path.exists(os.path.join(new_project_directory, 'pageGen.txt')):
        logger.info(
            'pageGen.txt already exists in %s ' % new_project_directory
        )
        return

    # Create ssh command.
    random_letters = ''.join([
        random.choice(string.ascii_lowercase + string.digits)
        for x in xrange(6)
     ])
    html_directory = os.path.join(
        config.get('start_makes', 'PageGenHtml'), new_name
    )
    ssh_command_string = 'ssh %s %s -p %s -u %s -s %s %s %s >> %s' % (
        config.get('start_makes', 'PageGenHost'),
        config.get('start_makes', 'PageGenPath'), id_number,
        investigator.lower(), random_letters + investigator.lower(),
        new_project_directory, html_directory,
        os.path.join(new_project_directory, 'pageGen.txt')
    )

    logger.info('ssh command string: %s' % ssh_command_string)
    # Extract username and password from standard output.
    pid = subprocess.Popen(
        ssh_command_string, shell=True, stdout=subprocess.PIPE
    )
    logger.info('ssh return value: %s' % pid.communicate()[0])


def email(recipient, content, config):
    """Send email with given content from EmailSender (in pathway.cfg) to
    recipient.
    """
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.ehlo()
    server.starttls()
    sender = config.get('Globals', 'EmailSender')
    msg = '\r\n'.join([
        'From: %s' % sender, 'To: %s' % recipient,
        'Subject: Script output', '', content
    ])
    server.login(sender, config.get('Globals', 'EmailPassword'))
    server.sendmail(sender, recipient, msg)
    server.quit()


def get_dir_info(project_directory_name):
    """Extract principal investigator and ID number from name of
    project directory. (Format is Project_[PI]_..._[ID].)
    """
    dir_name = os.path.basename(project_directory_name).split('_')
    return {'pi': dir_name[1], 'id': dir_name[-1]}


def get_new_file_name(project_directory, logger):
    """Return name of new file, which is specified in first line of
    newFileName.txt in project directory.
    """
    try:
        f = open(os.path.join(project_directory, 'newFileName.txt'))
    except IOError:
        logger.warn(
            'newFileName.txt not found in %s.' % project_directory
        )
        return None
    f_name = f.readline().strip()
    if not f_name:
        logger.warn(
            'newFileName.txt in %s is present, but empty.' % project_directory
        )
    return f_name


def make_links(original_directory, config, logger, seen_new_proj_dirs):
    """Create links. Run script at pageGenPath on pageGenHost."""
    logger.info('Setting up links for %s.' % original_directory)
    for d in glob.glob(original_directory + '/Project*'):
        logger.info('Setting up sublinks for %s.' % d)
        new_name = get_new_file_name(d, logger)
        if new_name is None:
            continue
        dir_info = get_dir_info(d)
        new_proj_dir = os.path.join(
            config.get('Globals', 'OutDirectory'), new_name
        )
        if not os.path.exists(new_proj_dir):
            logger.warn(
                'Cannot find new project directory %s.' % new_proj_dir
            )
        logger.info(
            'Setting up output html for project directory %s, PI %s, and ID '
            'number %s.' % (new_proj_dir, dir_info['pi'], dir_info['id'])
        )
        if new_proj_dir not in seen_new_proj_dirs:
            ssh_project(
                new_proj_dir, dir_info['pi'], dir_info['id'], new_name,
                config, logger
            )
            seen_new_proj_dirs.append(new_proj_dir)


def email_results(original_directory, config, logger):
    """Email results to EmailRecipient (in pathway.cfg). Return True if
    successful; False otherwise.
    """
    logger.info('Checking links for %s.' % original_directory)
    page_gen_is_present = True
    for proj_dir in glob.glob(original_directory + '/Project*'):
        logger.info('Checking sublinks for %s.' % proj_dir)
        new_name = get_new_file_name(proj_dir, logger)
        new_proj_dir = os.path.join(
            config.get('Globals', 'OutDirectory'), new_name
        )
        if not os.path.exists(new_proj_dir):
            logger.warn(
                'Cannot find new project directory %s.' % new_proj_dir
            )
        if os.path.exists(os.path.join(new_proj_dir, 'pageGen.txt')):
            with open(os.path.join(new_proj_dir, 'pageGen.txt')) as f:
                for email_content in f.readlines():
                    if email_content.startswith('INFO'):
                        logger.info(
                            'Emailing results. Email content: %s' %
                            email_content
                        )
                        email(
                            config.get('Globals', 'EmailRecipient'), email_content,
                            config
                        )
        else:
            logger.warn('pageGen.txt is missing from %s.' % new_proj_dir)
            page_gen_is_present = False
    return page_gen_is_present


def create_new_file_name(project_directory, plate_id):
    """Create file name and write it to newFileName.txt."""
    if os.path.exists(os.path.join(project_directory, 'newFileName.txt')):
        return
    dir_info = get_dir_info(project_directory)
    new_name_base = '%s-%s_%s_%s' % (
        datetime.utcnow().strftime('%Y%m%d'), plate_id, dir_info['pi'],
        dir_info['id']
    )
    with open(os.path.join(project_directory, 'newFileName.txt'), 'w') as f:
        f.write(new_name_base)


def rsync_files(original_directory, config, logger):
    """rsync project directories to OutDirectory (in pathway.cfg)."""
    proj_dirs = glob.glob(original_directory + '/Project*')
    plate_id = os.path.basename(
        os.path.dirname(original_directory)
    ).split('_')[2]
    for proj_dir in proj_dirs:
        create_new_file_name(proj_dir, plate_id)
    for proj_dir in proj_dirs:
        if not os.path.isdir(proj_dir):
            logger.warn('Cannot find project directory %s.' % proj_dir)
        logger.info('Running FASTQ Screen.')
        for sample_dir in glob.glob('%s/Sample*' % proj_dir):
            #logger.info('Running FASTQ Screen on sample %s.' % sample_dir)
            ret_code = subprocess.Popen(
                'perl /home/seqproc/fastq_screen_new/fastq_screen_new.pl --conf '
                '/home/seqproc/fastq_screen_new/fastq_screen.conf --outdir '
                '%s/fastq_screen_results %s --force' % (proj_dir, sample_dir),
                shell=True
            )
            time.sleep(120)
        time.sleep(1200)
        logger.info('FASTQ Screen completed.')
        new_name_base = get_new_file_name(proj_dir, logger)
        new_name = os.path.join(
            config.get('Globals', 'OutDirectory'), new_name_base
        )
        if not os.path.exists(new_name):
            os.makedirs(new_name)
        logger.info('Starting rsync from %s to %s.' % (proj_dir, new_name))
        ret_code = subprocess.Popen(
            'rsync -v -r -u %s %s' % (proj_dir, new_name), shell=True
        )
        logger.info(
            'rsync return code: %s' % ret_code.communicate()[0]
        )
        logger.info('Finished rsync from %s to %s.' % (proj_dir, new_name))
        new_proj_dir = os.path.join(new_name, os.path.basename(proj_dir)) # verify

        # add proj dir to to_be_aligned file
        to_be_aligned_fname = config.get('Globals', 'ToBeAlignedList')
        if os.path.exists(to_be_aligned_fname):
            with open(to_be_aligned_fname) as f:
                to_be_aligned_list = [z.strip() for z in f.readlines()]
        else:
            to_be_aligned_list = []
        if new_proj_dir not in to_be_aligned_list:
            to_be_aligned_list.append(new_proj_dir)
        with open(to_be_aligned_fname, 'w') as f:
            f.write('\n'.join(to_be_aligned_list))

        # Stone rsync
        dir_info = get_dir_info(new_proj_dir)
        if dir_info['pi'] == 'Stone' or (not int(dir_info['id'])):  # is a Stone lab run
            logger.info('Changing permissions on %s.' % new_proj_dir)
            subprocess.Popen(
                'find %s -type d -exec chmod 777 {} \;' % new_proj_dir, shell=True
            )
            subprocess.Popen(
                'find %s -type f -exec chmod 755 {} \;' % new_proj_dir, shell=True
            )
            new_name_stone = os.path.join(
                '/mnt/IIHG_transfers/', '%s_%s' %
                (new_name_base, os.path.basename(new_proj_dir))
            )
            if os.path.exists(new_name_stone):
                logger.info(
                    'Stone directory %s already exists. will not attempt '
                    'rsync.' % new_name_stone
                )
            else:
                logger.info(
                    'This is a Stone run rsync from %s to IVR %s.' %
                    (new_proj_dir, new_name_stone)
                )
                ret_code = subprocess.Popen(
                    'rsync -v -p -r %s %s' % (new_proj_dir, new_name_stone),
                    shell=True
                )
                logger.info(
                    'Stone rsync return code: %s.' % ret_code.communicate()[0]
                )
                email(
                    'adam-deluca@uiowa.edu',
                    'Finished rsync from %s to %s.' % (
                        new_proj_dir, new_name_stone
                    ),
                    config
                )
                
def has_required_files(dir_name, string_of_files, logger):
    """Return True if all files in comma-separated string_of_files are
    present in dir; otherwise False.
    """
    has_req_files = True
    for f in [z.strip() for z in string_of_files.split(',')]:
        f_path = os.path.join(dir_name, f)
        logger.info('Checking for %s.' % f_path)
        if not os.path.exists(f_path):
            has_req_files = False
    logger.info(
        '%s %s all required files.' %
        (dir_name, 'has' if has_req_files else 'does not have')
    )
    return has_req_files


def check_for_single_csv(dir_, logger):
    """
    If directory contains exactly one .csv file, return the path to the
    file; otherwise return None.
    """
    logger.info('Checking %s for sample sheets.' % dir_)
    csv_files = glob.glob('%s/*.csv' % dir_)
    if csv_files:
        logger.info('Found CSV files: %s.' % ','.join(csv_files))

    # Ignore processedSampleSheetN files.
    unprocessed_csv_files = [
        z for z in csv_files if not
        os.path.basename(z).startswith('processedSampleSheet')
    ]
    if len(unprocessed_csv_files) == 1:
        logger.info('Found single unprocessed CSV file in %s.' % dir_)
        return unprocessed_csv_files[0]
    else:
        if not unprocessed_csv_files:
            logger.info('Cannot find unprocessed sample sheet in %s.' % dir_)
        else:
            logger.info(
                'Found multiple unprocessed sample sheets in %s.' % dir_
            )
        return None


def has_correct_illumina_dir_form(dir_name, logger):
    """Return True if dir name has 3 underscores and ends with XX;
    otherwise False.
    """
    dir_name = os.path.basename(dir_name)
    has_correct_underscores = dir_name.count('_') == 3
    ends_with_XX = dir_name.endswith('XX')
    logger.info(
        'Directory %s: Number of underscores %s Illumina result directory '
        'format.' % (
            dir_name, 'satisfies' if has_correct_underscores else
            'does not satisfy'
        )
    )
    logger.info(
        'Directory %s: End of directory name %s Illlumina result directory '
        'format.' %
        (
            dir_name, 'satisfies' if ends_with_XX else
            'does not satisfy'
        )
    )
    looks_like_illumina_dir = has_correct_underscores and ends_with_XX
    logger.info(
        '%s %s an Illumina result directory.' % (
            dir_name, 'looks like' if looks_like_illumina_dir else
            'does not look like'
        )
    )
    return looks_like_illumina_dir


def process_sample_sheet(index_length, dir_name, input_csv, logger):
    """Given an input CSV, copy only the rows containing indices of given
    length into a new CSV file.
    """
    out_csv_path = os.path.join(
        dir_name, 'processedSampleSheet%s.csv' % index_length
    )
    if os.path.exists(out_csv_path):
        logger.warn(
            'Found already-existing processed sample sheet(s) at %s. '
            '(Directory was already run?)' % out_csv_path
        )
    in_csv_f = open(input_csv)
    with open(input_csv) as in_csv_f:
        with open(out_csv_path, 'w') as out_csv_f:
            first_line = in_csv_f.readline().strip('\n').strip('\r')
            if first_line.startswith('[Data]'):
                header = in_csv_f.readline().strip('\n').strip('\r')
            else:
                header = first_line
            # Assumes index col is called 'index' and
            # dual index second index col is called 'index2'
            index_col = header.split(',').index('index')
            index2_col = -1
            if 'index2' in header.split(','):
                index2_col = header.split(',').index('index2')
            out_csv_f.write('[Data]\n%s\n' % header)
            index_length_is_present = False
            for row in in_csv_f:
                if not row.strip():
                    continue
                # Change periods and spaces to _ so bcl2fastq doesn't choke.
                row = re.sub(r'[. ]', '_', row.strip())
                # Removed so that dual indexing works.
                # row = row.replace('-', '_')
                cols = row.strip('\n').strip('\r').split(',')
                curr_index_len = len(
                    [z for z in cols[index_col] if z not in '-_']
                )
                if index2_col != -1:
                    curr_index_len += len(
                        [z for z in cols[index2_col] if z not in '-_']
                    )
                if curr_index_len  == index_length:
                    index_length_is_present = True
                    out_csv_f.write('%s\n' % ','.join(cols))
    if index_length_is_present:
        logger.info(
            'Created processed sample sheet %s (index length %d).' %
            (out_csv_path, index_length)
        )
        return out_csv_path
    else:
        os.remove(out_csv_path)
        logger.info('No indices present of length %d.' % index_length)
        return None


def run_bcl2fastq(
    base_directory, index_length, sample_sheet, logger, dual_index_present
):
    """Run bcl2fastq for given index length. (Create Makefiles.)
    Return path to bcl2fastq output (Unaligned directory).
    """
    logger.info('Running bcl2fastq for index length %d.' % index_length)
    out_path = os.path.join(base_directory, 'Unaligned%d' % index_length)
    args = [
        #'/opt/illumina/bin/configureBclToFastq.pl',
        '/usr/local/bin/bcl2fastq',
        '--sample-sheet', sample_sheet,
        '--input-dir', '%s/Data/Intensities/BaseCalls' % base_directory,
        '--runfolder-dir', base_directory,
        '--output-dir', out_path,
        #'--ignore-missing-bcl',
        '--ignore-missing-bcls',
        '--ignore-missing-filter',
        '--ignore-missing-positions',
        '--demultiplexing-threads', '8',
        '--processing-threads', '16'
        #'--ignore-missing-stat',
        #'--fastq-cluster-count=0'
    ]
    if index_length in xrange(6, 9):
        if dual_index_present:
            bases_mask = 'y*,I%dN*,N*,y*' % index_length
        else:
            bases_mask = 'y*,I%dN*,y*' % index_length
        bcl_process = subprocess.Popen(
            args + ['--use-bases-mask=%s' % bases_mask],
            stdout=subprocess.PIPE
        )
    elif index_length > 8:
        bcl_process = subprocess.Popen(
            args + [
                '--use-bases-mask=y*,I%dN*,I%dN*,y*' %
                (index_length / 2, index_length / 2)
            ],
            stdout=subprocess.PIPE
        )
    else:
        bcl_process = subprocess.Popen(args, stdout=subprocess.PIPE)

    bcl_process.communicate()
    if not bcl_process.returncode:
        logger.info(
            'bcl2fastq completed successfully for index length %d.' %
            index_length
        )
    else:
        logger.info(
            'Error in bcl2fastq for index length %d. Return code: %d.' %
            (index_length, bcl_process.returncode)
        )

    return out_path


if __name__ == '__main__':
    sys.exit(main())

