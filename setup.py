from setuptools import setup,find_packages

setup(name='molSimplifyAD',version="v0.9-alpha",packages=find_packages(),
      entry_points={'console_scripts': ['msad = molSimplifyAD.__main__:main']},
      package_data={
          'molSimplifyAD':["molSimplifyAD/*.sh","molSimplifyAD/default_ligands.txt"]
      },
      include_package_data = True,
     )
