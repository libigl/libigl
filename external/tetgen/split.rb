#!/opt/local/bin/ruby

require 'fileutils'
FileUtils.mkdir_p 'src'

filename = 'tetgen.cxx';
out_filename = "header.txt"
out_file = File.open(out_filename, 'w') 

File.readlines(filename).each do |line|
  if line =~ %r!//// ([^_]*)_cxx!
    if out_filename != "src/#{$1}.cxx"
      out_filename = "src/#{$1}.cxx"
      puts out_filename
      out_file.close
      out_file = File.open(out_filename, 'w') 
      out_file.write("#include \"../tetgen.h\"\n")
    end
  end
  out_file.write(line)
end
out_file.close
