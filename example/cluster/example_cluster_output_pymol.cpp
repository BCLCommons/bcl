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
#include "cluster/bcl_cluster_output_pymol.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_dendrogram.h"
#include "cluster/bcl_cluster_input_table.h"
#include "cluster/bcl_cluster_linkage_average.h"
#include "cluster/bcl_cluster_node_colorer.h"
#include "cluster/bcl_cluster_node_description_average.h"
#include "math/bcl_math_const_function.h"
#include "util/bcl_util_color_gradient.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    // explicit instantiation
    template class NodeColorer< std::string, float>;
    template class NodeDescriptionAverage< std::string, float>;
    template class OutputPymolLabelString< float>;
    template class DistancesStored< std::string, float>;
  } // namespace cluster

  namespace math
  {
    // explicit instantiation
    template class ConstFunction< std::string, double>;
  } // namespace math

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_output_pymol.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterOutputPymol :
    public ExampleInterface
  {
  public:

    ExampleClusterOutputPymol *Clone() const
    {
      return new ExampleClusterOutputPymol( *this);
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
      util::BinaryFunctionSTLWrapper< std::less< float> >::s_Instance.IsDefined();
      const std::string input_table_filename( AddExampleInputPathToFilename( e_Cluster, "table_lower_triangle.txt"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, input_table_filename);
      cluster::InputTable< float> input_table( false, true);

      // create "member_distance_function" which will be used to calculate distances between members of a node
      // so that the center of the node can be calculated
      util::ShPtr
      <
        math::FunctionInterface
        <
          storage::VectorND< 2, util::SiPtr< const std::string> >, float
        >
      > member_distance_function( input_table.HandleInput( read));
      io::File::CloseClearFStream( read);

      util::ShPtr< util::BinaryFunctionInterface< float, float, bool> > compare( math::Comparisons< float>::GetEnums().e_Less);

      cluster::Dendrogram< std::string, float> dendrogram
      (
        util::ShPtr< cluster::LinkageInterface< std::string, float> >
        (
          new cluster::LinkageAverage< std::string, float>( member_distance_function)
        ),
        input_table.GetInputObjects()
      );

      util::ShPtr
      <
        util::FunctionInterface
        <
          storage::Triplet //!< argument type
          <
            util::SiPtr< const cluster::Node< std::string, float> >,
            storage::VectorND< 4, double>,
            std::string
          >,
          std::string //!< return type
        >
      > label_maker( new cluster::OutputPymolLabelString< float>( member_distance_function, compare, util::GetColors().e_Magenta));

      // color gradient function
      const util::ShPtr< util::FunctionInterface< double, linal::Vector3D> > color_gradient
      (
        new util::ColorGradient
        (
          math::Range< double>( 0, 10),
          storage::Vector< util::Color>::Create( util::GetColors().e_Red, util::GetColors().e_Yellow)
        )
      );
      // method for describing node. in this case returns a constant value of 5
      util::ShPtr< util::FunctionInterface< cluster::Node< std::string, float>, double> > node_describer
      (
        new cluster::NodeDescriptionAverage< std::string, float>
        (
          util::ShPtr< util::FunctionInterface< std::string, double> >
          (
            new math::ConstFunction< std::string, double>( 5.0)
          )
        )
      );
      util::ShPtr< util::FunctionInterface< cluster::Node< std::string, float>, linal::Vector3D> > node_colorer
      (
        new cluster::NodeColorer< std::string, float>( node_describer, color_gradient)
      );

      const std::string dendrogram_filename_correct
      (
        AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymol_correct.py")
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "default constructor");
      cluster::OutputPymol< std::string, float> output_format_protein_def;

      // parameter constructor
      BCL_MessageStd( "constructor taking parameters");

      cluster::OutputPymol< std::string, float> param_constr
      (
        100,
        50,
        100,
        node_colorer,
        250, //< number of labels to print
        label_maker,
        util::GetColors().e_Magenta,
        10,
        true,
        storage::VectorND< 2, float>( 0, 20),
        false
      );

      // check that the parameters were correctly set by the constructor
      BCL_ExampleCheck
      (
        param_constr.GetCylinderLength() == 100 &&
        param_constr.GetCylinderRadius() == 50 &&
        param_constr.GetCylinderSeparation() == 100 &&
        param_constr.GetNodeColorer() == node_colorer &&
        param_constr.GetMaxNodeLabels() == 250 &&
        param_constr.GetLabelMaker() == label_maker,
        true
      );

      // clone function
      BCL_MessageStd( "clone function");
      util::ShPtr< cluster::OutputPymol< std::string, float> > output_format_protein_clone
      (
        param_constr.Clone()
      );
      BCL_ExampleCheck( output_format_protein_clone.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        ( GetStaticClassName< cluster::OutputPymol< std::string, float> >()),
        output_format_protein_clone->GetClassIdentifier()
      );

      // SetCylinderLength changes "m_CylinderLength"
      output_format_protein_def.SetCylinderLength( 100);
      BCL_ExampleCheck( output_format_protein_def.GetCylinderLength(), 100);

      // SetCylinderRadius changes "m_CylinderRadius"
      output_format_protein_def.SetCylinderRadius( 50);
      BCL_ExampleCheck( output_format_protein_def.GetCylinderRadius(), 50);

      // SetCylinderSeparation changes "m_CylinderSeparation"
      output_format_protein_def.SetCylinderSeparation( 100);
      BCL_ExampleCheck( output_format_protein_def.GetCylinderSeparation(), 100);

      // GetNodeColorer changes "m_NodeColorer"
      output_format_protein_def.SetNodeColorer( node_colorer);
      BCL_ExampleCheck( output_format_protein_def.GetNodeColorer(), node_colorer);

      // SetMaxNodeLabels changes "m_MaxNodeLabels"
      output_format_protein_def.SetMaxNodeLabels( 250);
      BCL_ExampleCheck( output_format_protein_def.GetMaxNodeLabels(), 250);

      // SetLabelMaker changes "m_LabelMaker"
      output_format_protein_def.SetLabelMaker( label_maker);
      BCL_ExampleCheck( output_format_protein_def.GetLabelMaker(), label_maker);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test WriteOutput function
      {
        const std::string dendrogram_filename
        (
          AddExampleOutputPathToFilename( dendrogram, "dendrogram_OutputPymol.py")
        );

        param_constr.WriteOutput
        (
          dendrogram_filename,
          util::SiPtrList< const cluster::Node< std::string, float> >( 1, dendrogram.GetNode())
        );

        // check the dendrogram file is correct
        BCL_ExampleCheck( io::File::FilesMatch( dendrogram_filename, dendrogram_filename_correct), true);
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterOutputPymol

  const ExampleClass::EnumType ExampleClusterOutputPymol::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterOutputPymol())
  );

} // namespace bcl
