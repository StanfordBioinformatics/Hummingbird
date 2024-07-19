from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
        name='CloudHummingbird',
        version='1.1.0',
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
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
        ],
        python_requires='>=3.8, <4',
        install_requires=[
            'dsub==0.3.6',
            'future==0.18.3',
            'configparser==5.0.0',
            'scipy==1.10.0',
            'numpy==1.23.0',
            'matplotlib==3.5.2',
            'scikit-learn==1.5.0',
            'google-cloud-storage==1.30.0',
            'boto3==1.18.6',
            'azure-storage-blob==12.13.0',
            'azure-identity==1.16.1',
            'azure-batch==10.0.0',
            'azure-mgmt-compute==18.0.0',
            'retry==0.9.2',
        ],
        extras_require={
            'tests': ['mock']
        },
        include_package_data=True,
        packages=find_packages(
            exclude=('docs', 'scripts', 'Hummingbird/test')
        ),
        package_data={
            'Hummingbird': ['conf/examples/*'],
            'Hummingbird': ['AWS/*'],
            'Hummingbird': ['Azure/*'],
            'Hummingbird': ['*.md'],
            'Hummingbird': ['plot/*'],
            'Hummingbird': ['*.conf'],
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
