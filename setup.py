from setuptools import setup, find_packages

setup(
    name='kmersv',
    version='0.0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'kmersv=kmersv.kmersv:main'
        ]
    },
    install_requires=[
        'numpy>=1.20.1',
        'matplotlib>=3.7.0',
        'pandas>=1.2.0'
    ]
)
