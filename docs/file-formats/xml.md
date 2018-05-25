.xml - serialization format
===========================

------------------------------------------------------------------------

A .xml file contains the serialization of an object data structure generated with the XMLSerializer:

The top level elements represent the groups in which the object are organised. The object names are unique within these groups.

    <group1>
      <object1 val="value of object 1"/>
      <object2 val="value of object 2"/>
    </group1>

    <group2>
      <object1 val="value of object 1"/>
    </group2>

An object can be of following type:

- Basic types: **char, char\*, std::string, bool, usigned int, int, float, double**
- STL containers: **std::array, std::vector, std::pair**
- Eigen types: **Eigen::Matrix, Eigen::SparseMatrix**
- User defined types: **XMLSerializable\*.**

There can also be a hierarchical structure like `vector<int>`, this will result in the following serialization:

    <group>
      <vector size="3">
        <value0 val="1"/>
        <value1 val="2"/>
        <value2 val="3"/>
      </vector>
    </group>

An example of a serialization of an instance of the class Test

    class Test{
      int var1;
      vector&ltfloat> vec1;  
    };

is shown here:

    <group>
      <Test>
        <var1 val="0">
        <vec1 size="2">
          <value0 val="1"/>
          <value1 val="2"/>
        </vector>
      </Test>
    </group>

In the following we show the serialization of Eigen matrices.

`Eigen::Matrix<int,4,3>`:

    <group>
      <matrix row="4" col="3" matrix="
    1,2,3,
    4,5,6,
    7,8,9,
    10,11,12/>
    </group>

`Eigen::SparseMatrix<int>` (3x3 identity saved as triplets of the non-zero entries):

    <group>
      <matrix row="3" col="3" matrix="
    0,0,1,
    1,1,1,
    2,2,1/>
    </group>
