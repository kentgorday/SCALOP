from setuptools import setup
setup(name='scalop',
      version='1.0.0',
      description='Sequence-based antibody CAnonical LOoP structure annotation',
      author='Wing Ki Wong',
      author_email='opig@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk/webapps/SCALOP',
      packages=['scalop',
                'scalop.prosci',
                'scalop.prosci.util',
                'scalop.prosci.loops'],
      package_dir={'scalop': 'lib/python/scalop',
                   'scalop.prosci': 'lib/python/scalop/prosci',
                   'scalop.prosci.util': 'lib/python/scalop/prosci/util',
                   'scalop.prosci.loops': 'lib/python/scalop/prosci/loops'},
      package_data={
          'scalop': ['database/*'],
      },
      scripts=['bin/SCALOP'],
      license="BSD 3-clause"
     )
