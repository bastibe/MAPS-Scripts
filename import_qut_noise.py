import jbof
import pathlib
import numpy

import shutil
shutil.rmtree('QUT_NOISE', ignore_errors=True)

root = pathlib.Path('QUT-NOISE_original')

dataset = jbof.create_dataset('QUT_NOISE', {
    'README': (root / 'docs' / 'README.text').read_text(),
    'LICENSE': (root / 'QUT-NOISE' / 'LICENSE.txt').read_text()})

for noisefile in root.glob('QUT-NOISE/*.wav'):
    item = dataset.add_item(noisefile.stem, {})
    item.add_array_from_file('signal', noisefile, {})

    labelfile = (root / 'QUT-NOISE' / 'labels' / (noisefile.stem + '.lab.txt'))
    labels = []
    for line in labelfile.open():
        start, stop, label = line.split(maxsplit=2)
        labels.append((float(start), float(stop), label))
    labels = numpy.array(labels, dtype=[('start', float),
                                        ('stop', float),
                                        ('label', 'U32')])
    item.add_array('labels', labels, {})

    impulsefile = (root / 'QUT-NOISE' / 'impulses' / (noisefile.stem + '.imp.txt'))
    if impulsefile.exists():
        impulse = numpy.loadtxt(impulsefile)
        item.add_array('impulse', impulse, {})
