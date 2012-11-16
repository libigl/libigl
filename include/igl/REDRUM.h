#ifndef IGL_REDRUM_H
#define IGL_REDRUM_H

// Q: These should probably be inside the igl namespace. What's the correct
// way to do that?
// A: I guess the right way is to not use a macro but a proper function with
// streams as input and output.

// ANSI color codes for formating iostream style output

// Bold Red, etc.
#define REDRUM(X)      "\e[1m\e[31m"<<X<<"\e[m"
#define GREENRUM(X)    "\e[1m\e[32m"<<X<<"\e[m"
#define YELLOWRUM(X)   "\e[1m\e[33m"<<X<<"\e[m"
#define BLUERUM(X)     "\e[1m\e[34m"<<X<<"\e[m"
#define MAGENTARUM(X)  "\e[1m\e[35m"<<X<<"\e[m"
#define CYANRUM(X)     "\e[1m\e[36m"<<X<<"\e[m"

// Regular Red, etc.
#define REDGIN(X)      "\e[31m"<<X<<"\e[m"
#define GREENGIN(X)    "\e[32m"<<X<<"\e[m"
#define YELLOWGIN(X)   "\e[33m"<<X<<"\e[m"
#define BLUEGIN(X)     "\e[34m"<<X<<"\e[m"
#define MAGENTAGIN(X)  "\e[35m"<<X<<"\e[m"
#define CYANGIN(X)     "\e[36m"<<X<<"\e[m"

#endif 
