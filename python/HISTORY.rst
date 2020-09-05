=======
History
=======

0.1.0 (2019-05-21)
------------------

* First release on PyPI.

0.1.4 (2019-07-19)
------------------

* release planets as its own library

  * This means it's now a dependency of this project

* First documenation pages are now live at heat1d.readthedocs.org

0.3.1 (2020-09-04)
------------------

* rename classes to capital letters following best practices
* add Configurator object to be able to run `Model` with different chi values

  * This also allows change of gridding parameters and other settings
  * display a Configurator object like shown in the example to see its content
  
* depend on planets >= 0.4.6 due to bug in 0.4.5
* add plotting module and mpl stylesheet use