// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "cluster/bcl_cluster_output_pymol_label_string.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_dendrogram.h"
#include "cluster/bcl_cluster_input_table.h"
#include "cluster/bcl_cluster_linkage_average.h"
#include "cluster/bcl_cluster_node_colorer.h"
#include "cluster/bcl_cluster_output_pymol.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_output_pymol_label_string.cpp
  //!
  //! @author alexanns
  //! @date Sep 20, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace cluster
  {
    template class OutputPymolLabelString< double>;
  } // namespace cluster

  class ExampleClusterOutputPymolLabelString :
    public ExampleInterface
  {
  public:

    ExampleClusterOutputPymolLabelString *Clone() const
    {
      return new ExampleClusterOutputPymolLabelString( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      util::BinaryFunctionSTLWrapper< std::less< double> >::s_Instance.IsDefined();
      const std::string input_table_filename( AddExampleInputPathToFilename( e_Cluster, "table_lower_triangle.txt"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, input_table_filename);
      cluster::InputTable< double> input_table( false, true);

      // create "member_distance_function" which will be used to calculate distances between members of a node
      // so that the center of the node can be calculated
      util::ShPtr
      <
        math::FunctionInterface
        <
          storage::VectorND< 2, util::SiPtr< const std::string> >, double
        >
      > member_distance_function( input_table.HandleInput( read));
      io::File::CloseClearFStream( read);

      util::ShPtr< util::BinaryFunctionInterface< double, double, bool> > compare( math::Comparisons< double>::GetEnums().e_Less);

      cluster::Dendrogram< std::string, double> dendrogram
      (
        util::ShPtr< cluster::LinkageInterface< std::string, double> >
        (
          new cluster::LinkageAverage< std::string, double>( member_distance_function)
        ),
        input_table.GetInputObjects()
      );

      const std::string dendrogram_filename_correct
      (
        AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_correct.py")
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // call default constructor
      cluster::OutputPymolLabelString< double> def_constr;

      // parameter constructor
      const cluster::OutputPymolLabelString< double> param_constr
      (
        member_distance_function, compare, util::GetColors().e_Magenta, true
      );

      {
        const cluster::OutputPymol< std::string, double> output_pymol
        (
          200, 5, 10,
          util::ShPtr< util::FunctionInterface< cluster::Node< std::string, double>, linal::Vector3D> >( new cluster::NodeColorer< std::string, double>()),
          15,
          util::ShPtr
          <
            util::FunctionInterface
            <
              storage::Triplet //!< argument type
              <
                util::SiPtr< const cluster::Node< std::string, double> >,
                storage::VectorND< 4, double>,
                std::string
              >,
              std::string //!< return type
            >
          >( param_constr.Clone()),
          util::GetColors().e_Magenta,
          10,
          true,
          storage::VectorND< 2, double>( util::GetUndefinedDouble(), util::GetUndefinedDouble()),
          false
        );

        const std::string dendrogram_basename( AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_param_constr"));

        output_pymol.WriteOutput
        (
          dendrogram_basename + ".py",
          util::SiPtrList< const cluster::Node< std::string, double> >( 1, dendrogram.GetNode())
        );

        // check the dendrogram file is correct
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance
          (
            dendrogram_basename + ".py",
            dendrogram_filename_correct,
            0.001
          ),
          true
        );

        // check the list files are correct
        for( size_t list( 1); list <= 7; ++list)
        {
          if( list != 5)
          {
            BCL_ExampleCheck
            (
              io::File::FilesMatch
              (
                dendrogram_basename + "_node_" + util::Format()( list) + ".ls",
                AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_node_" + util::Format()( list) + "_correct.ls")
              ), true
            );
          }
        }
      }

      // test copy constructor
      {
        const cluster::OutputPymolLabelString< double> copy_constr( param_constr);

        const cluster::OutputPymol< std::string, double> output_pymol
        (
          200, 5, 10,
          util::ShPtr< util::FunctionInterface< cluster::Node< std::string, double>, linal::Vector3D> >( new cluster::NodeColorer< std::string, double>()),
          15,
          util::ShPtr
          <
            util::FunctionInterface
            <
              storage::Triplet //!< argument type
              <
                util::SiPtr< const cluster::Node< std::string, double> >,
                storage::VectorND< 4, double>,
                std::string
              >,
              std::string //!< return type
            >
          >( copy_constr.Clone()),
          util::GetColors().e_Magenta,
          10,
          true,
          storage::VectorND< 2, double>( util::GetUndefinedDouble(), util::GetUndefinedDouble()),
          false
        );

        const std::string dendrogram_basename( AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_copy_constr"));

        output_pymol.WriteOutput
        (
          dendrogram_basename + ".py",
          util::SiPtrList< const cluster::Node< std::string, double> >( 1, dendrogram.GetNode())
        );

        // check the dendrogram file is correct
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance
          (
            dendrogram_basename + ".py",
            dendrogram_filename_correct,
            0.001
          ),
          true
        );

        // check the list files are correct
        for( size_t list( 1); list <= 7; ++list)
        {
          if( list != 5)
          {
            BCL_ExampleCheck
            (
              io::File::FilesMatch
              (
                dendrogram_basename + "_node_" + util::Format()( list) + ".ls",
                AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_node_" + util::Format()( list) + "_correct.ls")
              ), true
            );
          }
        }
      }

      // virtual copy constructor
      {
        const util::ShPtr< cluster::OutputPymolLabelString< double> > clone_constr( param_constr.Clone());

        const cluster::OutputPymol< std::string, double> output_pymol
        (
          200, 5, 10,
          util::ShPtr< util::FunctionInterface< cluster::Node< std::string, double>, linal::Vector3D> >( new cluster::NodeColorer< std::string, double>()),
          15,
          clone_constr,
          util::GetColors().e_Magenta,
          10,
          true,
          storage::VectorND< 2, double>( util::GetUndefinedDouble(), util::GetUndefinedDouble()),
          false
        );

        const std::string dendrogram_basename( AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_clone_constr"));

        output_pymol.WriteOutput
        (
          dendrogram_basename + ".py",
          util::SiPtrList< const cluster::Node< std::string, double> >( 1, dendrogram.GetNode())
        );

        // check the dendrogram file is correct
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance
          (
            dendrogram_basename + ".py",
            dendrogram_filename_correct,
            0.001
          ),
          true
        );

        // check the list files are correct
        for( size_t list( 1); list <= 7; ++list)
        {
          if( list != 5)
          {
            BCL_ExampleCheck
            (
              io::File::FilesMatch
              (
                dendrogram_basename + "_node_" + util::Format()( list) + ".ls",
                AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelString_node_" + util::Format()( list) + "_correct.ls")
              ), true
            );
          }
        }
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        param_constr.GetClassIdentifier(),
        GetStaticClassName( cluster::OutputPymolLabelString< double>())
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // test operator()
      // operator() is tested above implicitly in the constructors to make sure it generates the correct text for the
      // dendrogram script and the node lists

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // test GetCylinderText
      BCL_ExampleCheck
      (
        cluster::OutputPymolLabelString< double>::GetCylinderText
        (
          "o_cluster_output_pymol_protein_model_0_0",
          0,
          3,
          57.1429,
          "6",
          0.05,
          0.4,
          util::GetColors().e_Magenta
        ),
        "axes=[[  0.0,0.0,0.4],[0.0,0.4,0.0],[0.0,0.0,0.0]]\n"
        "cyl_text(o_cluster_output_pymol_protein_model_0_0,plain,[0,3,57.1429],'6',0.05,color=Magenta,axes=axes)\n"
      );

      // test GetLabelCommand
      BCL_ExampleCheck
      (
        util::SplitString
        (
          cluster::OutputPymolLabelString< double>::GetLabelCommand
          (
            dendrogram.GetNode(),
            "test",
            0.0,
            10.0,
            5.0,
            2.0,
            util::GetColors().e_Magenta
          )
        ),
        util::SplitString
        (
          "test=[]\n"
          "axes=[[ 0.0,0.0,0.25],[0.0,0.25,0.0],[0.0,0.0,0.0]]\n"
          "cyl_text(test,plain,[0,6.25,1.20536],'4',0.03125,color=Magenta,axes=axes)\n"
          "axes=[[ 0.0,0.0,0.25],[0.0,0.25,0.0],[0.0,0.0,0.0]]\n"
          "cyl_text(test,plain,[0,5.9375,1.20536],'14',0.03125,color=Magenta,axes=axes)\n"
          "axes=[[ 0.0,0.0,0.25],[0.0,0.25,0.0],[0.0,0.0,0.0]]\n"
          "cyl_text(test,plain,[0,6.5625,1.20536],'5',0.03125,color=Magenta,axes=axes)\n"
          "labels_all.extend( test)"
        )
      );

      // test GetTextZCoordinate gives the z coordinate where the text for a label should start
      BCL_ExampleCheck( cluster::OutputPymolLabelString< double>::GetTextZCoordinate( 175.0, 50.0), 200.0);

      // test GetTextBlockSize gives the size that a piece of text should have
      BCL_ExampleCheck( cluster::OutputPymolLabelString< double>::GetTextBlockSize( 40.0, 20.0), 1.0);

      // test GetCylinderTextRadius gives the radius of the cylinders that create a CGO text
      BCL_ExampleCheck( cluster::OutputPymolLabelString< double>::GetCylinderTextRadius( 32.0), 4.0);

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterOutputPymolLabelProteinModel

  const ExampleClass::EnumType ExampleClusterOutputPymolLabelString::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterOutputPymolLabelString())
  );

} // namespace bcl
