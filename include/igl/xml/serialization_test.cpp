//
// Copyright (C) 2014 Christian Schüller <schuellchr@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//#ifndef IGL_SERIALIZATION_TEST_H
//#define IGL_SERIALIZATION_TEST_H

//#include <igl/Timer.h>
#include <igl/xml/serialize_xml.h>

namespace igl
{

  struct Test1111
  {
  };

  struct Test1 : public XMLSerializable
  {
    std::string ts;
    std::vector<Test1*> tvt;
    Test1* tt;

    Test1()
    {
      tt = NULL;
    }

    void InitSerialization()
    {
      Add(ts,"ts",false);
      Add(tvt,"tvt");
      Add(tt,"tt");
    }
  };

  struct Test2: public igl::XMLSerializableBase
  {
    char tc;
    int* ti;
    std::vector<short> tvb;
    float tf;

    Test2()
    {
      tc = '1';
      ti = NULL;
      tf = 1.0004;
      tvb.push_back(2);
      tvb.push_back(3);
    }

    void Serialize(tinyxml2::XMLDocument* doc,tinyxml2::XMLElement* element) const
    {
      igl::serialize_xml(tc,"tc",doc,element);
      igl::serialize_xml(ti,"ti",doc,element);
      igl::serialize_xml(tvb,"tvb",doc,element);
    }
    void Deserialize(const tinyxml2::XMLDocument* doc,const tinyxml2::XMLElement* element)
    {
      igl::deserialize_xml(tc,"tc",doc,element);
      igl::deserialize_xml(ti,"ti",doc,element);
      igl::deserialize_xml(tvb,"tvb",doc,element);
    }
    void Serialize(std::vector<char>& buffer) const
    {
      igl::serialize(tc,"tc",buffer);
      igl::serialize(ti,"ti",buffer);
      igl::serialize(tvb,"tvb",buffer);
      igl::serialize(tf,"tf",buffer);
    }
    void Deserialize(const std::vector<char>& buffer)
    {
      igl::deserialize(tc,"tc",buffer);
      igl::deserialize(ti,"ti",buffer);
      igl::deserialize(tvb,"tvb",buffer);
      igl::deserialize(tf,"tf",buffer);
    }
  };

  void serialization_test()
  {
    std::string file("test");

    bool tbIn = true,tbOut;
    char tcIn = 't',tcOut;
    unsigned char tucIn = 'u',tucOut;
    short tsIn = 6,tsOut;
    int tiIn = -10,tiOut;
    unsigned int tuiIn = 10,tuiOut;
    float tfIn = 1.0005,tfOut;
    double tdIn = 1.000000005,tdOut;

    int* tinpIn = NULL,*tinpOut = NULL;
    float* tfpIn = new float,*tfpOut = NULL;
    *tfpIn = 1.11101;

    std::string tstrIn("test12345"),tstrOut;

    Test2 tObjIn,tObjOut;
    int ti = 2;
    tObjIn.ti = &ti;


    Test1 test1,test2,test3;
    test1.ts = "100";
    test2.ts = "200";
    test3.ts = "300";

    Test1 testA, testC;
    testA.tt = &test1;
    testA.ts = "test123";
    testA.tvt.push_back(&test2);
    testA.tvt.push_back(&test3);

    Test1 testB = testA;
    testB.ts = "400";
    testB.tvt.pop_back();

    std::pair<int,bool> tPairIn(10,true);
    std::pair<int,bool> tPairOut;

    std::vector<int> tVector1In ={1,2,3,4,5};
    std::vector<int> tVector1Out;

    std::pair<int,bool> p1(10,1);
    std::pair<int,bool> p2(1,0);
    std::pair<int,bool> p3(10000,1);
    std::vector<std::pair<int,bool> > tVector2In ={p1,p2,p3};
    std::vector<std::pair<int,bool> > tVector2Out;

    std::set<std::pair<int,bool> > tSetIn ={p1,p2,p3};
    std::set<std::pair<int,bool> > tSetOut;

    std::map<int,bool> tMapIn ={p1,p2,p3};
    std::map<int,bool> tMapOut;

    Eigen::Matrix<float,3,3> tDenseMatrixIn;
    tDenseMatrixIn << Eigen::Matrix<float,3,3>::Random();
    tDenseMatrixIn.coeffRef(0,0) = 1.00001;
    Eigen::Matrix<float,3,3> tDenseMatrixOut;

    Eigen::Matrix<float,3,3,Eigen::RowMajor> tDenseRowMatrixIn;
    tDenseRowMatrixIn << Eigen::Matrix<float,3,3,Eigen::RowMajor>::Random();
    Eigen::Matrix<float,3,3,Eigen::RowMajor> tDenseRowMatrixOut;

    Eigen::SparseMatrix<double> tSparseMatrixIn;
    tSparseMatrixIn.resize(3,3);
    tSparseMatrixIn.insert(0,0) = 1.3;
    tSparseMatrixIn.insert(1,1) = 10.2;
    tSparseMatrixIn.insert(2,2) = 100.1;
    tSparseMatrixIn.finalize();
    Eigen::SparseMatrix<double> tSparseMatrixOut;

    // binary serialization

    igl::serialize(tbIn,file);
    igl::deserialize(tbOut,file);
    assert(tbIn == tbOut);

    igl::serialize(tcIn,file);
    igl::deserialize(tcOut,file);
    assert(tcIn == tcOut);

    igl::serialize(tucIn,file);
    igl::deserialize(tucOut,file);
    assert(tucIn == tucOut);

    igl::serialize(tsIn,file);
    igl::deserialize(tsOut,file);
    assert(tsIn == tsOut);

    igl::serialize(tiIn,file);
    igl::deserialize(tiOut,file);
    assert(tiIn == tiOut);

    igl::serialize(tuiIn,file);
    igl::deserialize(tuiOut,file);
    assert(tuiIn == tuiOut);

    igl::serialize(tfIn,file);
    igl::deserialize(tfOut,file);
    assert(tfIn == tfOut);

    igl::serialize(tdIn,file);
    igl::deserialize(tdOut,file);
    assert(tdIn == tdOut);

    igl::serialize(tinpIn,file);
    igl::deserialize(tinpOut,file);
    assert(tinpIn == tinpOut);

    igl::serialize(tfpIn,file);
    igl::deserialize(tfpOut,file);
    assert(*tfpIn == *tfpOut);
    tfpOut = NULL;

    igl::serialize(tstrIn,file);
    igl::deserialize(tstrOut,file);
    assert(tstrIn == tstrOut);

    // updating
    igl::serialize(tbIn,"tb",file,true);
    igl::serialize(tcIn,"tc",file);
    igl::serialize(tiIn,"ti",file);
    tiIn++;
    igl::serialize(tiIn,"ti",file);
    tiIn++;
    igl::serialize(tiIn,"ti",file);
    igl::deserialize(tbOut,"tb",file);
    igl::deserialize(tcOut,"tc",file);
    igl::deserialize(tiOut,"ti",file);
    assert(tbIn == tbOut);
    assert(tcIn == tcOut);
    assert(tiIn == tiOut);

    igl::serialize(tsIn,"tsIn",file,true);
    igl::serialize(tVector1In,"tVector1In",file);
    igl::serialize(tVector2In,"tsIn",file);
    igl::deserialize(tVector2Out,"tsIn",file);
    for(unsigned int i=0;i<tVector2In.size();i++)
    {
      assert(tVector2In[i].first == tVector2Out[i].first);
      assert(tVector2In[i].second == tVector2Out[i].second);
    }
    tVector2Out.clear();

    igl::serialize(tObjIn,file);
    igl::deserialize(tObjOut,file);
    assert(tObjIn.tc == tObjOut.tc);
    assert(*tObjIn.ti == *tObjOut.ti);
    for(unsigned int i=0;i<tObjIn.tvb.size();i++)
      assert(tObjIn.tvb[i] == tObjOut.tvb[i]);
    tObjOut.ti = NULL;

    igl::serialize(tPairIn,file);
    igl::deserialize(tPairOut,file);
    assert(tPairIn.first == tPairOut.first);
    assert(tPairIn.second == tPairOut.second);

    igl::serialize(tVector1In,file);
    igl::deserialize(tVector1Out,file);
    for(unsigned int i=0;i<tVector1In.size();i++)
      assert(tVector1In[i] == tVector1Out[i]);

    igl::serialize(tVector2In,file);
    igl::deserialize(tVector2Out,file);
    for(unsigned int i=0;i<tVector2In.size();i++)
    {
      assert(tVector2In[i].first == tVector2Out[i].first);
      assert(tVector2In[i].second == tVector2Out[i].second);
    }

    igl::serialize(tSetIn,file);
    igl::deserialize(tSetOut,file);
    assert(tSetIn.size() == tSetOut.size());

    igl::serialize(tMapIn,file);
    igl::deserialize(tMapOut,file);
    assert(tMapIn.size() == tMapOut.size());

    igl::serialize(tDenseMatrixIn,file);
    igl::deserialize(tDenseMatrixOut,file);
    assert((tDenseMatrixIn - tDenseMatrixOut).sum() == 0);

    igl::serialize(tDenseRowMatrixIn,file);
    igl::deserialize(tDenseRowMatrixOut,file);
    assert((tDenseRowMatrixIn - tDenseRowMatrixOut).sum() == 0);

    igl::serialize(tSparseMatrixIn,file);
    igl::deserialize(tSparseMatrixOut,file);
    assert((tSparseMatrixIn - tSparseMatrixOut).sum() == 0);

    igl::serialize(testB,file);
    igl::deserialize(testC,file);
    assert(testB.ts == testC.ts);
    assert(testB.tvt.size() == testC.tvt.size());
    for(unsigned int i=0;i<testB.tvt.size();i++)
    {
      assert(testB.tvt[i]->ts == testC.tvt[i]->ts);
      assert(testB.tvt[i]->tvt.size() == testC.tvt[i]->tvt.size());
      assert(testB.tvt[i]->tt == testC.tvt[i]->tt);
    }
    assert(testB.tt->ts == testC.tt->ts);
    assert(testB.tt->tvt.size() == testC.tt->tvt.size());
    assert(testB.tt->tt == testC.tt->tt);
    testC = Test1();

    // big data test
    /*std::vector<std::vector<float> > bigDataIn,bigDataOut;
    for(unsigned int i=0;i<10000;i++)
    {
    std::vector<float> v;
    for(unsigned int j=0;j<10000;j++)
    {
    v.push_back(j);
    }
    bigDataIn.push_back(v);
    }

    igl::Timer timer;
    timer.start();
    igl::serialize(bigDataIn,file);
    timer.stop();
    std::cout << "ser: " << timer.getElapsedTimeInMilliSec() << std::endl;

    timer.start();
    igl::deserialize(bigDataOut,file);
    timer.stop();
    std::cout << "des: " << timer.getElapsedTimeInMilliSec() << std::endl;
    char c;
    std::cin >> c; */

    // xml serialization

    igl::serialize_xml(tbIn,file);
    igl::deserialize_xml(tbOut,file);
    assert(tbIn == tbOut);

    igl::serialize_xml(tcIn,file);
    igl::deserialize_xml(tcOut,file);
    assert(tcIn == tcOut);

    igl::serialize_xml(tucIn,file);
    igl::deserialize_xml(tucOut,file);
    assert(tucIn == tucOut);

    igl::serialize_xml(tsIn,file);
    igl::deserialize_xml(tsOut,file);
    assert(tsIn == tsOut);

    igl::serialize_xml(tiIn,file);
    igl::deserialize_xml(tiOut,file);
    assert(tiIn == tiOut);

    igl::serialize_xml(tuiIn,file);
    igl::deserialize_xml(tuiOut,file);
    assert(tuiIn == tuiOut);

    igl::serialize_xml(tfIn,file);
    igl::deserialize_xml(tfOut,file);
    assert(tfIn == tfOut);

    igl::serialize_xml(tdIn,file);
    igl::deserialize_xml(tdOut,file);
    assert(tdIn == tdOut);

    igl::serialize_xml(tinpIn,file);
    igl::deserialize_xml(tinpOut,file);
    assert(tinpIn == tinpOut);

    igl::serialize_xml(tfpIn,file);
    igl::deserialize_xml(tfpOut,file);
    assert(*tfpIn == *tfpOut);

    igl::serialize_xml(tstrIn,file);
    igl::deserialize_xml(tstrOut,file);
    assert(tstrIn == tstrOut);

    // updating
    igl::serialize_xml(tbIn,"tb",file,false,true);
    igl::serialize_xml(tcIn,"tc",file);
    igl::serialize_xml(tiIn,"ti",file);
    tiIn++;
    igl::serialize_xml(tiIn,"ti",file);
    tiIn++;
    igl::serialize_xml(tiIn,"ti",file);
    igl::deserialize_xml(tbOut,"tb",file);
    igl::deserialize_xml(tcOut,"tc",file);
    igl::deserialize_xml(tiOut,"ti",file);
    assert(tbIn == tbOut);
    assert(tcIn == tcOut);
    assert(tiIn == tiOut);

    igl::serialize_xml(tsIn,"tsIn",file,false,true);
    igl::serialize_xml(tVector1In,"tVector1In",file);
    igl::serialize_xml(tVector2In,"tsIn",file);
    igl::deserialize_xml(tVector2Out,"tsIn",file);
    for(unsigned int i=0;i<tVector2In.size();i++)
    {
      assert(tVector2In[i].first == tVector2Out[i].first);
      assert(tVector2In[i].second == tVector2Out[i].second);
    }
    tVector2Out.clear();

    // binarization
    igl::serialize_xml(tVector2In,"tVector2In",file,true);
    igl::deserialize_xml(tVector2Out,"tVector2In",file);
    for(unsigned int i=0;i<tVector2In.size();i++)
    {
      assert(tVector2In[i].first == tVector2Out[i].first);
      assert(tVector2In[i].second == tVector2Out[i].second);
    }

    igl::serialize_xml(tObjIn,file);
    igl::deserialize_xml(tObjOut,file);
    assert(tObjIn.tc == tObjOut.tc);
    assert(*tObjIn.ti == *tObjOut.ti);
    for(unsigned int i=0;i<tObjIn.tvb.size();i++)
      assert(tObjIn.tvb[i] == tObjOut.tvb[i]);

    igl::serialize_xml(tPairIn,file);
    igl::deserialize_xml(tPairOut,file);
    assert(tPairIn.first == tPairOut.first);
    assert(tPairIn.second == tPairOut.second);

    igl::serialize_xml(tVector1In,file);
    igl::deserialize_xml(tVector1Out,file);
    for(unsigned int i=0;i<tVector1In.size();i++)
      assert(tVector1In[i] == tVector1Out[i]);

    igl::serialize_xml(tVector2In,file);
    igl::deserialize_xml(tVector2Out,file);
    for(unsigned int i=0;i<tVector2In.size();i++)
    {
      assert(tVector2In[i].first == tVector2Out[i].first);
      assert(tVector2In[i].second == tVector2Out[i].second);
    }

    igl::serialize_xml(tSetIn,file);
    igl::deserialize_xml(tSetOut,file);
    assert(tSetIn.size() == tSetOut.size());

    igl::serialize_xml(tMapIn,file);
    igl::deserialize_xml(tMapOut,file);
    assert(tMapIn.size() == tMapOut.size());

    igl::serialize_xml(tDenseMatrixIn,file);
    igl::deserialize_xml(tDenseMatrixOut,file);
    assert((tDenseMatrixIn - tDenseMatrixOut).sum() == 0);

    igl::serialize_xml(tDenseRowMatrixIn,file);
    igl::deserialize_xml(tDenseRowMatrixOut,file);
    assert((tDenseRowMatrixIn - tDenseRowMatrixOut).sum() == 0);

    igl::serialize_xml(tSparseMatrixIn,file);
    igl::deserialize_xml(tSparseMatrixOut,file);
    assert((tSparseMatrixIn - tSparseMatrixOut).sum() == 0);

    igl::serialize_xml(testB,file);
    igl::deserialize_xml(testC,file);
    assert(testB.ts == testC.ts);
    assert(testB.tvt.size() == testC.tvt.size());
    for(unsigned int i=0;i<testB.tvt.size();i++)
    {
      assert(testB.tvt[i]->ts == testC.tvt[i]->ts);
      assert(testB.tvt[i]->tvt.size() == testC.tvt[i]->tvt.size());
      assert(testB.tvt[i]->tt == testC.tvt[i]->tt);
    }
    assert(testB.tt->ts == testC.tt->ts);
    assert(testB.tt->tvt.size() == testC.tt->tvt.size());
    assert(testB.tt->tt == testC.tt->tt);

    // big data test
    /*std::vector<std::vector<float> > bigDataIn,bigDataOut;
    for(unsigned int i=0;i<10000;i++)
    {
    std::vector<float> v;
    for(unsigned int j=0;j<10000;j++)
    {
    v.push_back(j);
    }
    bigDataIn.push_back(v);
    }

    igl::Timer timer;
    timer.start();
    igl::serialize_xml(bigDataIn,"bigDataIn",file,igl::SERIALIZE_BINARY);
    timer.stop();
    std::cout << "ser: " << timer.getElapsedTimeInMilliSec() << std::endl;

    timer.start();
    igl::deserialize_xml(bigDataOut,"bigDataIn",file);
    timer.stop();
    std::cout << "des: " << timer.getElapsedTimeInMilliSec() << std::endl;
    char c;
    std::cin >> c;*/

    std::cout << "All tests run successfully!\n";
  }
}

//#endif
