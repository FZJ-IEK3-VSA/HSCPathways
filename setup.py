from setuptools import setup, find_packages

setup(
    name='HSCPathways',
    version='3.1.4',
    description='Evaluation of hydrogen infrastructure.',
 
    url='http://www.fz-juelich.de/iek/iek-3/DE/Home/home_node.html',
    author='Markus Reuss',
    author_email='m.reuss@fz-juelich.de',    
    license='',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'CoolProp',
        'jupyter',
        'openpyxl',
        'xlrd',
        'XlsxWriter',
        'pandas',
        'matplotlib==1.5.1',
    ]
)