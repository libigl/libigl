#include <test_common.h>
#include <igl/dirname.h>

#include <string>
#include <vector>
#include <tuple>

TEST_CASE("dirname: examples", "[igl]")
{
  const std::vector<
    std::tuple<std::string,std::string> >
    examples = 
  {
    std::make_tuple("/foo"                 ,"/"           ),
    std::make_tuple("/foo/"                ,"/"           ),
    std::make_tuple("/foo//"               ,"/"           ),
    std::make_tuple("/foo/./"              ,"/foo"        ),
    std::make_tuple("/foo/bar"             ,"/foo"        ),
    std::make_tuple("/foo/bar."            ,"/foo"        ),
    std::make_tuple("/foo/bar.txt"         ,"/foo"        ),
    std::make_tuple("/foo/bar.txt.zip"     ,"/foo"        ),
    std::make_tuple("/foo/bar.dir/"        ,"/foo"        ),
    std::make_tuple("/foo/bar.dir/file"    ,"/foo/bar.dir"),
    std::make_tuple("/foo/bar.dir/file.txt","/foo/bar.dir"),
    std::make_tuple("."                   ,"."           ),
    std::make_tuple("../foo"              ,".."           ),
    std::make_tuple("./foo"              ,"."           ),
    std::make_tuple("../foo/"             ,".."           ),
    std::make_tuple("foo/"               ,"."           ),
    std::make_tuple("foo//"               ,"."           ),
    std::make_tuple("foo/./"              ,"foo"        ),
    std::make_tuple("foo/bar"             ,"foo"        ),
    std::make_tuple("foo/bar."            ,"foo"        ),
    std::make_tuple("foo/bar.txt"         ,"foo"        ),
    std::make_tuple("foo/bar.txt.zip"     ,"foo"        ),
    std::make_tuple("foo/bar.dir/"        ,"foo"        ),
    std::make_tuple("foo/bar.dir/file"    ,"foo/bar.dir"),
    std::make_tuple("foo/bar.dir/file.txt","foo/bar.dir")
  };
  for(const auto & example : examples)
  {
    std::string d;
    d = igl::dirname(std::get<0>(example));
    REQUIRE (d == std::get<1>(example));
  }
}

