from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

_DEPENDENCIES = [
        'dsub==0.3.6',
        'future==0.18.2',
        'configparser==5.0.0',
        'scipy==1.5.2',
        'numpy==1.19.1',
        'matplotlib==3.3.0',
        'scikit-learn==0.23.2',
        'google-cloud-storage==1.30.0',
        'boto3==1.14.38',
]
setup(
        name='CloudHummingbird',
        version='1.0.1',
        description='A tool that recommends the best way to run your genomics pipelines on the cloud',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/StanfordBioinformatics/Hummingbird',
        author='Stanford Center For Genomics And Personalized Medicine',
        author_email='abahman@stanford.edu',
        classifiers=[
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 2.7',
        ],
        package_dir={'': 'src'},
        packages=find_packages(where='src'),
        python_requires='>=2.7, <4',
        install_requires = _DEPENDENCIES,
        #extras_require={
        #    'aws': ['boto3==1.14.38'],
        #    'gcp': ['google-cloud-storage==1.30.0','dsub==0.3.6'],
        #},
        include_package_data=True,
        package_data={
            'Hummingbird':['conf/examples/*'],
            'Hummingbird':['AWS/*'],
            'Hummingbird':['*.md'],
            'Hummingbird':['plot/*'],
            'Hummingbird':['*.conf'],
        },
        entry_points={
            'console_scripts': [
                'hummingbird=Hummingbird.hummingbird:main',
            ],
        },
        project_urls={
            'Bug Reports': 'https://github.com/StanfordBioinformatics/Hummingbird/issues',
            'Source': 'https://github.com/StanfordBioinformatics/Hummingbird',
        },
)
