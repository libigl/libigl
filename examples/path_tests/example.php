<?php
  if($argc <= 1)
  {
    printf("USAGE:\n  php example.php [path_1] [path_2] ... [path_n]\n");
    return 1;
  }
  for($i = 1; $i < $argc; $i++)
  {
    printf("file_exists(%s) --> %s\n",$argv[$i],(file_exists($argv[$i])?"TRUE":"FALSE"));
    printf("is_dir(%s) --> %s\n",$argv[$i],(is_dir($argv[$i])?"TRUE":"FALSE"));
    printf("is_file(%s) --> %s\n",$argv[$i],(is_file($argv[$i])?"TRUE":"FALSE"));
    printf("is_readable(%s) --> %s\n",$argv[$i],(is_readable($argv[$i])?"TRUE":"FALSE"));
    printf("is_writable(%s) --> %s\n",$argv[$i],(is_writable($argv[$i])?"TRUE":"FALSE"));
  }
?>
