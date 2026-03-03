"""
conftest.py - Automatically add the build directory to sys.path so that
the 'geometry' extension module can be imported.

The module is located in <repo_root>/build/python/ after building with CMake.
The build directory is searched relative to the repository root.
"""

import sys
import os

# Walk up from the tests directory to find the repo root
_tests_dir = os.path.dirname(__file__)
_repo_root = os.path.abspath(os.path.join(_tests_dir, '..', '..'))

# Support overriding the build directory via GEOMETRY_BUILD_DIR environment variable
_build_dir = os.environ.get(
    'GEOMETRY_BUILD_DIR',
    os.path.join(_repo_root, 'build', 'python')
)

if _build_dir not in sys.path:
    sys.path.insert(0, _build_dir)
