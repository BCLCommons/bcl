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
#include "cluster/bcl_cluster_output_pymol_label_small_molecule.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_dendrogram.h"
#include "cluster/bcl_cluster_input_table.h"
#include "cluster/bcl_cluster_linkage_classes.h"
#include "cluster/bcl_cluster_node_colorer.h"
#include "cluster/bcl_cluster_output_pymol.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_output_pymol_label_small_molecule.cpp
  //! @details this example test the ClusterOutputPymolLabelSmallMolecule class
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace cluster
  {
    template class OutputPymolLabelSmallMolecule< double>;
  } // namespace cluster

  class ExampleClusterOutputPymolLabelSmallMolecule :
    public ExampleInterface
  {
  public:

    ExampleClusterOutputPymolLabelSmallMolecule *Clone() const
    {
      return new ExampleClusterOutputPymolLabelSmallMolecule( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      const std::string input_table_filename( AddExampleInputPathToFilename( e_Cluster, "small_mol_diff.table"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, input_table_filename);
      cluster::InputTable< double> input_table( true, false);

      // create "member_distance_function" which will be used to calculate distances between members of a node
      // so that the center of the node can be calculated
      util::ShPtr< cluster::DistancesStored< std::string, double> > member_distance_function( input_table.HandleInput( read));

      util::ShPtr< util::BinaryFunctionInterface< double, double, bool> > compare( math::Comparisons< double>::GetEnums().e_Less);

      BCL_MessageDbg( "number of input objects is " + util::Format()( input_table.GetInputObjects()->GetSize()));

      cluster::Dendrogram< std::string, double> dendrogram
      (
        util::ShPtr< cluster::LinkageInterface< std::string, double> >( new cluster::LinkageAverage< std::string, double>( member_distance_function)),
        input_table.GetInputObjects()
      );

      // the file which contains the correct dendrogram against which the other outputted dendrograms will be compared
      const std::string dendrogram_filename_correct
      (
        AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelSmallMolecule_correct.py")
      );

      // mdl filenames
      std::string mdl_files[] =
      {
        AddExampleOutputPathToFilename( dendrogram, "0.mol"),
        AddExampleOutputPathToFilename( dendrogram, "1.mol"),
        AddExampleOutputPathToFilename( dendrogram, "2.mol"),
        AddExampleOutputPathToFilename( dendrogram, "3.mol"),
        AddExampleOutputPathToFilename( dendrogram, "4.mol")
      };

      // node files
      std::string node_files[] =
      {
        AddExampleOutputPathToFilename( dendrogram, "node_1_0.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_2_1.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_3_2.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_4_3.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_5_4.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_7_1.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_7_4.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_8_0.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_8_3.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_9_0.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_9_1.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_9_3.mol"),
        AddExampleOutputPathToFilename( dendrogram, "node_9_4.mol")
      };

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const cluster::OutputPymolLabelSmallMolecule< double> def_constr;

      // parameter constructor
      cluster::OutputPymolLabelSmallMolecule< double> param_constr
      (
        member_distance_function,
        compare,
        AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"),
        util::GetColors().e_Magenta,
        true
      );
      {
        const cluster::OutputPymol< std::string, double> output_pymol
        (
          200, 5, 10,
          util::ShPtr< util::FunctionInterface< cluster::Node< std::string, double>, linal::Vector3D> >
            ( new cluster::NodeColorer< std::string, double>()),
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

        const std::string dendrogram_basename
        (
          AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelSmallMolecule_param_constr")
        );

        output_pymol.WriteOutput
        (
          dendrogram_basename + ".py",
          util::SiPtrList< const cluster::Node< std::string, double> >( 1, dendrogram.GetNode())
        );

        // check the dendrogram file is correct
        BCL_ExampleCheck( io::File::FilesMatch( dendrogram_basename + ".py", dendrogram_filename_correct), true);

        // check the mdl files are correct
        for( size_t file( 0); file < 5; ++file)
        {
          BCL_ExampleCheck( io::File::FilesMatch( mdl_files[ file], mdl_files[ file] + "_correct"), true);
        }

        // check the node member files are correct
        for( size_t file( 0); file < 13; ++file)
        {
          BCL_ExampleCheck( io::File::FilesMatch( node_files[ file], node_files[ file] + "_correct"), true);
        }
      }

      // test clone constructor
      {
        util::ShPtr< cluster::OutputPymolLabelSmallMolecule< double> > clone_constr( param_constr.Clone());
        const cluster::OutputPymol< std::string, double> output_pymol
        (
          200, 5, 10,
          util::ShPtr< util::FunctionInterface< cluster::Node< std::string, double>, linal::Vector3D> >
            ( new cluster::NodeColorer< std::string, double>()),
          15,
          clone_constr,
          util::GetColors().e_Magenta,
          10,
          true,
          storage::VectorND< 2, double>( util::GetUndefinedDouble(), util::GetUndefinedDouble()),
          false
        );

        const std::string dendrogram_basename
        (
          AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymolLabelSmallMolecule_clone_constr")
        );

        output_pymol.WriteOutput
        (
          dendrogram_basename + ".py",
          util::SiPtrList< const cluster::Node< std::string, double> >( 1, dendrogram.GetNode())
        );

        // check the dendrogram file is correct
        BCL_ExampleCheck( io::File::FilesMatch( dendrogram_basename + ".py", dendrogram_filename_correct), true);

        // check the mdl files are correct
        for( size_t file( 0); file < 5; ++file)
        {
          BCL_ExampleCheck( io::File::FilesMatch( mdl_files[ file], mdl_files[ file] + "_correct"), true);
        }

        // check the node member files are correct
        for( size_t file( 0); file < 13; ++file)
        {
          BCL_ExampleCheck( io::File::FilesMatch( node_files[ file], node_files[ file] + "_correct"), true);
        }
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( param_constr.GetClassIdentifier(), GetStaticClassName( param_constr));

    ///////////////
    // operators //
    ///////////////

      // test operator() is checked above in the constructors

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // test LoadPymolObject
      {
        // create const string "given_text" and initialize with the text given by the function
        const storage::Vector< std::string> given_text
        (
          util::SplitString
          (
            cluster::OutputPymolLabelSmallMolecule< double>::LoadPymolObject( "test", "test_obj")
          )
        );

        // create const string "correct_text" and initialize with the text that should be given by the function
        const storage::Vector< std::string> correct_text
        (
          util::SplitString( std::string( "cmd.load( \"test\", \"test_obj\")\n"))
        );

        // make sure that "correct_text" is the same as "given_text"
        BCL_ExampleCheck( given_text, correct_text);
      }

      // test TranslatePymolObject
      {
        // create const string "given_text" and initialize with the text given by the function
        const storage::Vector< std::string> given_text
        (
          util::SplitString
          (
            cluster::OutputPymolLabelSmallMolecule< double>::TranslatePymolObject
            (
              0.0, 1.0, 2.0, "test_obj"
            )
          )
        );

        // create const string "correct_text" and initialize with the text that should be given by the function
        const storage::Vector< std::string> correct_text
        (
          util::SplitString( std::string( "cmd.translate( [0, 1, 2], \"test_obj\")\n"))
        );

        // make sure that "correct_text" is the same as "given_text"
        BCL_ExampleCheck( given_text, correct_text);
      }

      // test GetObjectNameLabel
      {
        // create const string "given_text" and initialize with the text given by the function
        const storage::Vector< std::string> given_text
        (
          util::SplitString
          (
            cluster::OutputPymolLabelSmallMolecule< double>::GetObjectNameLabel
            (
              0.0, 2.0, 1.0, dendrogram.GetNode().GetCenter( member_distance_function, compare), 5.0, util::GetColors().e_Magenta
            )
          )
        );

        // create const string "correct_text" and initialize with the text that should be given by the function
        const storage::Vector< std::string> correct_text
        (
          util::SplitString
          (
            std::string( "center_0=[]\n") +
            std::string( "axes=[[ 0.0,0.0,0.05],[0.0,0.05,0.0],[0.0,0.0,0.0]]\n") +
            std::string( "cyl_text(center_0,plain,[0,1.75,2.86964],'0',0.00625,color=Magenta,axes=axes)\n") +
            std::string( "labels_all.extend(\n") +
            std::string( "center_0)\n")
          )
        );

        // make sure that "correct_text" is the same as "given_text"
        BCL_ExampleCheck( given_text, correct_text);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterOutputPymolLabelSmallMolecule

  const ExampleClass::EnumType ExampleClusterOutputPymolLabelSmallMolecule::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterOutputPymolLabelSmallMolecule())
  );

} // namespace bcl
