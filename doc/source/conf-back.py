import sys

sys.path.append( "/fs/home/mr937/breathe-master/breathe/" )

extensions = ['sphinx.ext.pngmath', 'sphinx.ext.todo', 'breathe' ]

breathe_projects = { "SEING": "/fs/home/mr937/SEING/doc/xml/" }

breathe_default_project = "SEING"