cd NL
cat nl_linkage.h nl.h |egrep -v "NL/" >> ../nl.h
echo "#include \"nl.h\"" > ../nl_single_file.c
cat `cat filelist.txt` | egrep -v "NL/" >> ../nl_single_file.c
cd ..


# for each source plugin specified, add the sources
#cd plugins
#    for plugin_dir in $*
#    do
#        cd ${plugin_dir}
#            echo Adding plugin ${plugin_dir}
#            #TODO
#        cd ..
#    done  
