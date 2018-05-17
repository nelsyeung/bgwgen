from setuptools import setup

setup(name='bgwgen',
      version='0.1.0',
      description='Generate BerkeleyGW simulation files from a single input',
      url='http://github.com/nelsyeung/bgwgen',
      author='Nelson Yeung',
      author_email='nelsyeung@icloud.com',
      license='MIT',
      packages=['bgwgen'],
      zip_safe=False,
      scripts=['bgwgen/bgwgen'],
      classifiers=[
          'Intended Audience :: Developers',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ])
