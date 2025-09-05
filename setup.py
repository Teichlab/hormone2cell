from setuptools import setup, find_packages

setup(
    name='hormone2cell',
    version='0.1.0',
    description='Predict hormone producing/receiving strength in single cell datasets',
    url='https://github.com/Teichlab/hormone2cell',
    packages=find_packages(exclude=['docs', 'notebooks']),
    install_requires=[
        'scanpy',
        'matplotlib-venn'
    ],
    package_data={
        "hormone2cell": ["*.pkl","*.h5ad"]
    },
    author='Lijiang Fei, Krzysztof Polanskiu',
    author_email='lf529@cam.ac.uk',
    license='non-commercial license'
)
