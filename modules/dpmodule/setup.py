from distutils.core import setup, Extension

def main():
    setup(name="dp_module",
          version="1.0.0",
          description="C implementation for DP in SVanchor.",
          author="Eru Shinonome",
          author_email="snnm.eru@gmail.com",
          ext_modules=[Extension("dp_module", ["dpmodule.c"])])

if __name__ == "__main__":
    main()
