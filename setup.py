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
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8'
        ],
        python_requires='>=3.6, <4',
        install_requires=[
            'dsub==0.3.6',
            'future==0.18.2',
            'configparser==5.0.0',
            'scipy==1.5.2',
            'numpy==1.19.1',
            'matplotlib==3.3.0',
            'scikit-learn==0.23.2',
            'google-cloud-storage==1.30.0',
            'boto3==1.18.6',
            'azure-storage-blob==12.6.0',
            'azure-identity==1.5.0',
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
