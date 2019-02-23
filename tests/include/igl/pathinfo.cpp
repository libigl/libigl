#include <test_common.h>
#include <igl/pathinfo.h>

#include <string>
#include <vector>
#include <tuple>

TEST_CASE("pathinfo: examples", "[igl]")
{
  const std::vector<
    std::tuple<std::string,std::string,std::string,std::string,std::string> >
    examples = 
  {
    std::make_tuple("/foo"                 ,"/"           ,"foo"        ,""   ,"foo"),
    std::make_tuple("/foo/"                ,"/"           ,"foo"        ,""   ,"foo"),
    std::make_tuple("/foo//"               ,"/"           ,"foo"        ,""   ,"foo"),
    std::make_tuple("/foo/./"              ,"/foo"        ,"."          ,""   ,""),
    std::make_tuple("/foo/bar"             ,"/foo"        ,"bar"        ,""   ,"bar"),
    std::make_tuple("/foo/bar."            ,"/foo"        ,"bar."       ,""   ,"bar"),
    std::make_tuple("/foo/bar.txt"         ,"/foo"        ,"bar.txt"    ,"txt","bar"),
    std::make_tuple("/foo/bar.txt.zip"     ,"/foo"        ,"bar.txt.zip","zip","bar.txt"),
    std::make_tuple("/foo/bar.dir/"        ,"/foo"        ,"bar.dir"    ,"dir","bar"),
    std::make_tuple("/foo/bar.dir/file"    ,"/foo/bar.dir","file"       ,""   ,"file"),
    std::make_tuple("/foo/bar.dir/file.txt","/foo/bar.dir","file.txt"   ,"txt","file")
  };
  for(const auto & example : examples)
  {
    std::string d,b,e,f;
    igl::pathinfo(std::get<0>(example),d,b,e,f);
    REQUIRE (d == std::get<1>(example));
    REQUIRE (b == std::get<2>(example));
    REQUIRE (e == std::get<3>(example));
    REQUIRE (f == std::get<4>(example));
  }
}
