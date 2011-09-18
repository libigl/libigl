<?php
  $file_handle = fopen("input.txt", "rb");
  while (!feof($file_handle) )
  {
    $string = trim(fgets($file_handle));
    if(!feof($file_handle) )
    {
      $p = pathinfo($string);
      echo "$string -> ".
        $p['dirname'].",".
        $p['basename'].",".
        $p['extension'].",".
        $p['filename']."\n";
    }
  }
  fclose($file_handle);
?>
