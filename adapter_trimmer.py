import re
import os
import PySimpleGUI as sg
import threading
import time
import queue

sg.theme('DarkTeal9')  # Change the theme here

# Layout

def create_layout():
    layout = [
        [
            sg.Text('Adapter file:', size=(15, 1), tooltip='Choose a FASTQ or FASTA file with adapter sequences'),
            sg.Input(key='adapter_file', tooltip='Choose a FASTQ or FASTA file with adapter sequences'),
            sg.FileBrowse(file_types=(('FASTQ Files', '*.fastq;*.fq'), ('FASTA Files', '*.fasta;*.fa')))
        ],
        [
            sg.Text('Sequence file:', size=(15, 1), tooltip='Choose a FASTQ file with sequences to be trimmed'),
            sg.Input(key='sequence_file', tooltip='Choose a FASTQ file with sequences to be trimmed'),
            sg.FileBrowse(file_types=(('FASTQ Files', '*.fastq;*.fq'),))
        ],
        [
            sg.Text('Output file:', size=(15, 1), tooltip='Choose the output FASTQ file to save trimmed sequences'),
            sg.Input(key='output_file', tooltip='Choose the output FASTQ file to save trimmed sequences'),
            sg.SaveAs(file_types=(('FASTQ Files', '*.fastq;*.fq'),))
        ],
        [
            sg.Button('Start Trimming', tooltip='Trim adapters from sequences'),
            sg.Button('Clear', tooltip='Clear input fields'),
            sg.Button('Help', tooltip='Show usage instructions'),
            sg.Button('Exit', tooltip='Exit the program')
        ],
        [sg.ProgressBar(100, orientation='h', size=(40, 20), key='progress_bar')],
        [sg.Text('', key='progress_text', size=(30, 1))],
        [sg.Output(size=(80, 20))],
        [sg.Column([
            [sg.Text('Adapter Trimmer', key='Application name', size=(None, 1), justification='left', font=("Alike", 11, "bold"))],  # Application name
            [sg.Text('Thesis Project. Created by Mohit Panwar. Supervised by Julia Åkesson.', key='credits', size=(None, 1), justification='left', font=("Alike", 9))]  # Credits
        ])]
    ]
    return layout

# Read adapter sequences from a file

def read_adapter_sequences(adapter_file):
    if not os.path.exists(adapter_file):
        raise ValueError(f'Adapter file {adapter_file} does not exist')
    if not os.access(adapter_file, os.R_OK):
        raise ValueError(f'Adapter file {adapter_file} is not readable')

    adapter_list = []
    _, adapter_format = os.path.splitext(adapter_file)
    adapter_format = adapter_format.lower().lstrip('.')
    with open(adapter_file) as f:
        if adapter_format in ['fasta', 'fa']:
            for count, line in enumerate(f, start=0):
                if count % 2 == 1:
                    adapter_list.append(line.strip())
        elif adapter_format in ['fastq', 'fq']:
            for count, line in enumerate(f, start=0):
                if count % 4 == 1:
                    adapter_list.append(line.strip())
        else:
            raise ValueError('Unrecognized file extension')

    return adapter_list

# Trim adapters from sequences

def trim_adapters(queue, error_queue, adapter_list, sequence_file, output_file, progress_callback=None):
    try:
        total_sequences = sum(1 for _ in open(sequence_file)) // 4
        processed_sequences = 0
        trimmed_sequences = 0
        start_time = time.time()

        with open(sequence_file) as seq_file, open(output_file, 'w') as out_file:
            patterns = [re.compile(adapter) for adapter in adapter_list]
            for count, line in enumerate(seq_file, start=0):
                if count % 4 == 0:
                    heading = line.strip()
                elif count % 4 == 1:
                    sequence = line.strip()
                elif count % 4 == 2:
                    plus_line = line.strip()
                elif count % 4 == 3:
                    quality_line = line.strip()

                    original_sequence = sequence
                    for pattern in patterns:
                        sequence = pattern.sub('', sequence)

                    if original_sequence != sequence:
                        trimmed_sequences += 1

                    out_file.write(f'{heading}\n{sequence}\n{plus_line}\n{quality_line}\n')

                    if progress_callback and count % 4 == 3:
                        processed_sequences += 1
                        progress = (processed_sequences / total_sequences) * 100
                        progress_callback(progress)

        elapsed_time = time.time() - start_time
        return trimmed_sequences, elapsed_time
    
    except Exception as e:
        error_queue.put(str(e))



# Main function

def main():
    window = sg.Window('Adapter Trimmer', create_layout())

    def update_progress_bar(progress):
        window['progress_bar'].update(progress)
        window['progress_text'].update(f'{int(progress * total_sequences / 100)} / {total_sequences} sequences processed')

    trimming_thread = None
    result_queue = queue.Queue()
    error_queue = queue.Queue()

    while True:
        event, values = window.read(timeout=100)

        if event == sg.WINDOW_CLOSED or event == 'Exit':
            break

        if event == 'Start Trimming':
            adapter_file = values['adapter_file']
            sequence_file = values['sequence_file']
            output_file = values['output_file']

            if not adapter_file:
                sg.popup('Please choose an adapter file')
                continue

            if not sequence_file:
                sg.popup('Please choose a sequence file')
                continue

            if not output_file:
                sg.popup('Please choose an output file')
                continue

            try:
                adapter_list = read_adapter_sequences(adapter_file)
                with open(output_file, 'w') as out_f:
                    out_f.write('')
                
                total_sequences = sum(1 for _ in open(values['sequence_file'])) // 4
                trimming_thread = threading.Thread(target=lambda q, eq, *args: q.put(trim_adapters(*args)), args=(result_queue, error_queue, adapter_list, sequence_file, output_file, update_progress_bar))
                trimming_thread.start()
            except Exception as e:
                sg.popup(f'Error: {e}')

        if event == 'Clear':
            window['adapter_file']('')
            window['sequence_file']('')
            window['output_file']('')

        if event == 'Help':
            help_text = """How to use Adapter Trimmer:
1. Choose an adapter file (FASTQ or FASTA format) containing the adapter sequences to be trimmed.
2. Choose a sequence file (FASTQ format) containing the sequences to be trimmed.
3. Choose an output file (FASTQ format) where the trimmed sequences will be saved.
4. Click "Start Trimming" to start the trimming process. A progress bar will indicate the progress of the operation.
5. When trimming is complete, a confirmation message will be displayed.

Note: You can click "Clear" to reset the input fields and start over."""
            sg.popup('Help', help_text)

        if trimming_thread and not trimming_thread.is_alive():
            trimmed_sequences, elapsed_time = result_queue.get()
            print(f'\nTrimming complete.\nTrimmed sequences: {trimmed_sequences}\nRuntime: {elapsed_time:.2f} seconds')
            trimming_thread = None

        if not error_queue.empty():
            error = error_queue.get()
            sg.popup(f'Trimming Error: {error}')
            error_queue.queue.clear()
            trimming_thread = None

    window.close()

if __name__ == '__main__':
    main()