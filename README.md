# A relaxation sim

This project can be built with [fpm](https://fpm.fortran-lang.org/index.html) using the following command:
```
fpm build relaxation --flag="-Ofast -march=native -mtune=native -flto -ftree-vectorize" 
```

The location of the executable can be found using:
```
find ./ -name "relaxation" | grep app
```

The plots can be viewed using the provided python file by changing the file-names in the list.
