Generates a figurecaption each Image which stands alone in a paragraph,
similar to [pandoc’s handling of images/figures](http://pandoc.org/README.html#extension-implicit_figures)

--------------------------------------------

Licensed under the GPL 2 (see LICENSE.md)

Copyright 2015 - Jan Dittrich by
building upon the [markdown-figures](https://github.com/helderco/markdown-figures) Plugin by
Copyright 2013 - [Helder Correia](http://heldercorreia.com) (GPL2)

--------------------------------------------

Example – this source:

    Bla bla bla

    ![this is the caption](http://lorempixel.com/400/200/)

    Next paragraph starts here

would generate this:

    <p> Bla bla bla</p>
    
    <figure>
        <img src="http://lorempixel.com/400/200/">
        <figcaption>this is the caption</figcaption>
    </figure>
   
    <p>Next paragraph starts here</p>
