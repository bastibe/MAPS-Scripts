import jbof
import pathlib
from collections import defaultdict
import re
import numpy

import shutil
shutil.rmtree('PTDB_TUG', ignore_errors=True)

root = pathlib.Path('PTDB-TUG_original')

dataset = jbof.create_dataset('PTDB_TUG', {
    'Recording Protocol': (root / 'RECORDING-PROTOCOL.txt').read_text(),
    'Speaker Profiles': (root / 'SPEAKER-PROFILES.txt').read_text(),
    'TIMIT Prompts': (root / 'TIMIT-PROMPTS.txt').read_text()})

sentences = {}
pattern = re.compile('([^.?!]+[.?!]).*\(([a-z0-9]+)\)')
for line in dataset.metadata['TIMIT Prompts'].split('\n'):
    if line.startswith(';') or not line:
        continue
    line, label = pattern.match(line).groups()
    sentences[label] = line

speaker_profiles = defaultdict(dict)
pattern = re.compile('([MF][0-9]{2})\s+([0-9]{2})\s+(Male|Female)\s+'
                     '(Ireland|USA|Canada|England|South Africa)\s+'
                     '(sa1,2 sx[0-9]+-[0-9]+\s+si[0-9]+-[0-9]+)\s*'
                     '(.*)')
for line in dataset.metadata['Speaker Profiles'].split('\n'):
    if line.startswith('Speaker') or line.startswith('-') or not line:
        continue
    speaker, age, sex, country, sentence, comment = pattern.match(line).groups()
    speaker_profiles[speaker] = {'Age': age, 'Sex': sex, 'Home Country': country,
                                 'Sentences': sentence, 'Comment': comment}


for path in root.glob('**/*.wav'):
    kind, speaker, sentence = path.stem.split('_')
    kind = dict(lar='laryngograph', mic='signal')[kind]
    print('\r', kind, speaker, sentence, end='')
    itemname = f'{speaker}_{sentence}'
    if not dataset.has_item(itemname):
        item = dataset.add_item(itemname, dict(speaker_profiles[speaker],
                                               speaker_id=speaker,
                                               sentence=sentences[sentence],
                                               sentence_id=sentence))
    else:
        item = dataset.get_item(itemname)
    item.add_array_from_file(kind, path, {})

for path in root.glob('**/*.f0'):
    kind, speaker, sentence = path.stem.split('_')
    print('\r', kind, speaker, sentence, end='')
    itemname = f'{speaker}_{sentence}'
    item = dataset.get_item(itemname)
    pitch = numpy.loadtxt(path)[:, 0]
    duration = len(item.signal) / item.signal.metadata['samplerate']
    sample_interval = duration / len(pitch)
    time = numpy.arange(len(pitch)) * sample_interval
    item.add_array('pitch', numpy.rec.fromarrays([time, pitch], names=['time', 'pitch']),
                   {'samplerate': 1/sample_interval})

print('') # newline
