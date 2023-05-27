import os
import sys
import PySimpleGUI as sg
import threading
import time
import queue

sg.theme('DarkTeal9')  # Change the theme here

def create_layout():
    layout = [
        [sg.Text('Sequence file:', size=(15, 1)), sg.Input(key='-SEQUENCE_FILE-'), sg.FileBrowse(file_types=(('FASTQ Files', '*.fastq;*.fq'),))],
        [sg.Text('Quality score threshold:', size=(15, 1)), sg.Input(key='-THRESHOLD-', size=(5,1))],
        [sg.Text('Output file:', size=(15, 1)), sg.Input(key='-OUTPUT_FILE-'), sg.FileSaveAs(file_types=(('FASTQ Files', '*.fastq;*.fq'),))],
        [sg.Button('Start Trimming'), sg.Button('Clear'), sg.Button('Help'), sg.Button('Exit')],
        [sg.ProgressBar(100, orientation='h', size=(40, 20), key='progress_bar')],
        [sg.Text('', key='progress_text', size=(30, 1))],
        [sg.Output(size=(80, 20))],
        [sg.Text('', key='result_text')],
        [sg.Text('Quality Trimmer', key='Application name', size=(None, 1), justification='left', font=("Alike", 11, "bold"))],  # Application name
        [sg.Text('Thesis Project. Created by Mohit Panwar. Supervised by Julia Åkesson.', key='credits', size=(None, 1), justification='left', font=("Alike", 9))]  # Credits
    ]
    return layout


def quality_trimmer(result_queue, error_queue, sequence_file, threshold, output_file, total_count_queue, progress_queue):
    try:
        if not os.path.exists(sequence_file):
            raise ValueError(f'Sequence file {sequence_file} does not exist')
        if not os.access(sequence_file, os.R_OK):
            raise ValueError(f'Sequence file {sequence_file} is not readable')

        base_file_name = os.path.splitext(os.path.basename(sequence_file))[0]
        if output_file is None:
            output_file = 'quality_trimmed_' + base_file_name
        else:
            output_dir, output_name = os.path.split(output_file)
            if not os.path.exists(output_dir):
                raise ValueError(f'Output directory {output_dir} does not exist')
            if os.path.isdir(output_file):
                output_file = os.path.join(output_file, 'quality_trimmed_' + base_file_name)
            elif os.path.splitext(output_name)[1] not in ['.fastq', '.fq']:
                output_file = output_file + '.fastq'

        trimmed_count = 0
        total_count = 0
        processed_count = 0
        start_time = time.time()

        with open(sequence_file, 'r') as f, open(output_file, 'w') as g:
            for line in f:
                total_count += 1
                if total_count % 4 == 1:
                    identifier = line.strip()
                elif total_count % 4 == 2:
                    sequence = line.strip()
                elif total_count % 4 == 3:
                    separator = line.strip()
                else:
                    quality_scores = [ord(c) - 33 for c in line.strip()]
                    trim_point = len(quality_scores)
                    for i in reversed(range(len(quality_scores))):
                        if quality_scores[i] >= threshold:
                            trim_point = i
                            break
                    trimmed_sequence = sequence[:trim_point+1]
                    trimmed_quality_scores = quality_scores[:trim_point+1]
                    if len(trimmed_sequence) > 0:
                        trimmed_count += 1
                        g.write(identifier + '\n')
                        g.write(trimmed_sequence + '\n')
                        g.write(separator + '\n')
                        g.write(''.join([chr(q + 33) for q in trimmed_quality_scores]) + '\n')
                    processed_count += 1
                    progress = processed_count / (total_count // 4) * 100
                    progress_queue.put(progress)

        total_count //= 4
        total_count_queue.put(total_count)

        discarded_count = total_count - trimmed_count
        discarded_percent = discarded_count / total_count * 100
        elapsed_time = time.time() - start_time
        result_queue.put((threshold, total_count, trimmed_count, discarded_count, discarded_percent, elapsed_time, output_file))
    except Exception as e:
        error_queue.put(str(e))

def main():
    window = sg.Window('Quality Trimmer', create_layout())
    progress_queue = queue.Queue()
    total_count_queue = queue.Queue()
    total_count = 0
    trimming_thread = None
    result_queue = queue.Queue()
    error_queue = queue.Queue()

    def update_progress_bar(progress):
        window['progress_bar'].update(progress)
        window['progress_text'].update(f'{int(progress * total_count / 100)} / {total_count} sequences processed')

    while True:
        event, values = window.read(timeout=100)

        if event == sg.WINDOW_CLOSED or event == 'Exit':
            break

        if event == 'Start Trimming':
            sequence_file = values['-SEQUENCE_FILE-']
            threshold = values['-THRESHOLD-']
            output_file = values['-OUTPUT_FILE-']

            if not sequence_file:
                sg.popup('Please choose a sequence file')
                continue

            if not threshold:
                sg.popup('Please enter a quality score threshold')
                continue

            if not output_file:
                sg.popup('Please choose an output file')
                continue

            try:
                threshold = int(threshold)
                with open(output_file, 'w') as out_f:
                    out_f.write('')
                
                trimming_thread = threading.Thread(target=quality_trimmer, args=(result_queue, error_queue, sequence_file, threshold, output_file, total_count_queue, progress_queue))
                trimming_thread.start()
            except Exception as e:
                sg.popup(f'Error: {e}')

        if trimming_thread and not trimming_thread.is_alive():
            result = result_queue.get()
            print(f'Trimming complete.\nTrimmed sequences: {result[2]}\nRuntime: {result[5]:.2f} seconds')
            trimming_thread = None

        if event == 'Clear':
            window['-SEQUENCE_FILE-'].update('')
            window['-THRESHOLD-'].update('')
            window['-OUTPUT_FILE-'].update('')

        if event == 'Help':
            help_text = """How to use Quality Trimmer:
1. Choose a sequence file (FASTQ format) containing the sequences to be trimmed.
2. Enter the quality score threshold for trimming.
3. Choose an output file (FASTQ format) where the trimmed sequences will be saved.
4. Click "Start Trimming" to start the trimming process. A progress bar will indicate the progress of the operation.
5. When trimming is complete, a confirmation message will be displayed.

Note: You can click "Clear" to reset the input fields and start over."""
            sg.popup('Help', help_text)

        if not total_count_queue.empty():
            total_count = total_count_queue.get()

        if not error_queue.empty():
            error = error_queue.get()
            sg.popup(f'Trimming Error: {error}')
            error_queue.queue.clear()
            trimming_thread = None

        if not progress_queue.empty():
            progress = progress_queue.get()
            update_progress_bar(progress)

    window.close()

if __name__ == '__main__':
    main()
