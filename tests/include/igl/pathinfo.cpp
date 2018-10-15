#include <test_common.h>
#include <igl/pathinfo.h>

#include <string>
#include <vector>
#include <tuple>

TEST(pathinfo, examples)
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
    ASSERT_EQ(std::get<1>(example),d);
    ASSERT_EQ(std::get<2>(example),b);
    ASSERT_EQ(std::get<3>(example),e);
    ASSERT_EQ(std::get<4>(example),f);
  }
}
