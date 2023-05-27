import os
import threading
import time
import queue
import PySimpleGUI as sg

sg.theme('DarkTeal9')

def create_layout():
    layout = [
        [sg.Text('Sequence file:', size=(15, 1)), sg.Input(tooltip="Select a '.fastq' or '.fq' file", key='-SEQUENCE_FILE-'), sg.FileBrowse(file_types=(('FASTQ Files', '*.fastq;*.fq'),))],
        [sg.Text('Quality score threshold:', size=(20, 1)), sg.Input(tooltip="Enter the quality score threshold", key='-THRESHOLD-', size=(5,1))],
        [sg.Text('Output file:', size=(15, 1)), sg.Input(tooltip="Specify the output file location", key='-OUTPUT_FILE-'), sg.FileSaveAs(file_types=(('FASTQ Files', '*.fastq;*.fq'),))],
        [sg.Button('Start Filtering'), sg.Button('Clear'), sg.Button('Help'), sg.Button('Exit')],
        [sg.Output(size=(80, 20))],
        [sg.Text('', key='result_text')],
        [sg.Text('Quality Filter', key='Application name', size=(None, 1), justification='left', font=("Alike", 11, "bold"))],
        [sg.Text('Thesis Project. Created by Mohit Panwar. Supervised by Julia Åkesson.', key='credits', size=(None, 1), justification='left', font=("Alike", 9))]
    ]
    return layout

def count_sequences(sequence_file):
    with open(sequence_file) as f:
        return sum(1 for line in f if line.startswith('@'))  # Count the number of lines starting with '@'

def quality_filter(sequence_file, threshold, output_file, total_sequences, progress_queue):
    filtered_count = 0
    total_count = 0
    start_time = time.time()
    with open(sequence_file, 'r') as f, open(output_file, 'w') as g:
        for line in f:
            identifier = line.strip()
            sequence = next(f).strip()
            separator = next(f).strip()
            quality_scores = [ord(c) - 33 for c in next(f).strip()]

            # We count a sequence every time we encounter an identifier line (which starts with '@')
            if identifier.startswith('@'):
                total_count += 1

            if sum(quality_scores)/len(quality_scores) < threshold:
                continue

            filtered_count += 1
            g.write(identifier + '\n')
            g.write(sequence + '\n')
            g.write(separator + '\n')
            g.write(''.join([chr(q + 33) for q in quality_scores]) + '\n')

    elapsed_time = time.time() - start_time
    discarded_count = total_count - filtered_count
    discarded_percent = discarded_count / total_count * 100
    progress_queue.put_nowait(('Result', (threshold, total_count, filtered_count, discarded_count, discarded_percent, elapsed_time, output_file)))




def main():
    window = sg.Window('Quality Filter', create_layout())
    total_count = 0
    progress_queue = queue.Queue()
    filtering_thread = None

    while True:
        event, values = window.read(timeout=100)

        if event == sg.WINDOW_CLOSED or event == 'Exit':
            break

        if event == 'Start Filtering':
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

            if not os.path.exists(sequence_file):
                sg.popup('Sequence file does not exist')
                continue
            if not os.access(sequence_file, os.R_OK):
                sg.popup('Sequence file is not readable')
                continue

            base_file_name = os.path.splitext(os.path.basename(sequence_file))[0]
            if output_file is None:
                output_file = 'quality_filtered_' + base_file_name
            else:
                output_dir, output_name = os.path.split(output_file)
                if not os.path.exists(output_dir):
                    sg.popup('Output directory does not exist')
                    continue
                if os.path.isdir(output_file):
                    output_file = os.path.join(output_file, 'quality_filtered_' + base_file_name)
                elif os.path.splitext(output_name)[1] not in ['.fastq', '.fq']:
                    output_file = output_file + '.fastq'

            try:
                threshold = int(threshold)
                total_count = count_sequences(sequence_file)
                filtering_thread = threading.Thread(target=quality_filter, args=(sequence_file, threshold, output_file, total_count, progress_queue), daemon=True)
                filtering_thread.start()
            except ValueError:
                sg.popup('Error: Threshold must be an integer.')
            except Exception as e:
                sg.popup(f'Error during filtering start: {e}')

        elif event == 'Clear':
            window['-SEQUENCE_FILE-'].update('')
            window['-THRESHOLD-'].update('')
            window['-OUTPUT_FILE-'].update('')
            window['Output'].update('')
            window['result_text'].update('')

        elif event == 'Help':
            sg.popup("This tool filters low-quality reads from a '.fastq' or '.fq' file based on the provided quality score threshold.\n\n1. Select a FASTQ file.\n2. Set a quality score threshold.\n3. Specify an output file.\n4. Click 'Start Filtering' to start the process.\n\nResults will be displayed in the output window after filtering is complete.")

        try:
            msg_type, msg_data = progress_queue.get_nowait()
            if msg_type == 'Result':
                threshold, total_count, filtered_count, discarded_count, discarded_percent, elapsed_time, output_file = msg_data
                print(f'Filtering completed in {elapsed_time:.2f} seconds.\nQuality cut-off: {threshold}\nInput: {total_count} reads\nOutput: {filtered_count} reads\nDiscarded: {discarded_count} reads ({discarded_percent:.2f}%)\nFiltered file is saved as: {output_file}')
            elif msg_type == 'Error':
                print('Error during quality filtering:', msg_data)
        except queue.Empty:
            pass

    window.close()

if __name__ == '__main__':
    main()
