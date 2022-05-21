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

#ifndef BCL_APP_CLUSTER_H_
#define BCL_APP_CLUSTER_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "cluster/bcl_cluster.fwd.hh"
#include "command/bcl_command.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Cluster
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_app_cluster.cpp @endlink
    //! @author alexanns
    //! @date Sep 17, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Cluster :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! input file with distances between objects
      util::ShPtr< command::FlagInterface> m_InputFile;

      //! format of input file
      util::ShPtr< command::FlagInterface> m_InputFormat;

      //! output format
      util::ShPtr< command::FlagInterface> m_OutputFormat;

      //! height cutoff
      util::ShPtr< command::FlagInterface> m_HeightCutoff;

      //! the output file for the results of the clustering
      util::ShPtr< command::FlagInterface> m_OutputFile;

      //! the linkage type to be used for clustering
      util::ShPtr< command::FlagInterface> m_Linkage;

      //! the distance type to be used for clustering - "<" for dissimilarity measures; ">" for similarity measures
      util::ShPtr< command::FlagInterface> m_DistanceType;

      //! output the dendrogram in a Python script that can be displayed by Pymol
      util::ShPtr< command::FlagStatic> m_OutputPymol;
      //! the unit length of the cylinders used to build the dendrogram
      util::ShPtr< command::ParameterInterface> m_CylinderLength;
      //! the radius of the cylinders used to build the dendrogram
      util::ShPtr< command::ParameterInterface> m_CylinderRadius;
      //! the unit separation between cylinders i.e. the base leaves will be separated by this amount
      util::ShPtr< command::ParameterInterface> m_CylinderSeparation;
      //! the maximum number of labels that should be displayed on the dendrogram
      util::ShPtr< command::ParameterInterface> m_MaxNodeLabels;
      //! the number of labels that should be along the vertical axis
      util::ShPtr< command::ParameterInterface> m_NumberLabelsAlongYScale;
      //! the name of the file which will have the Python code for displaying the dendrogram in Pymol
      util::ShPtr< command::ParameterInterface> m_OutputPymolOutputFile;

      //! specifies that nodes should be colored according to some description (e.g.rmsd to native; restraint agreement)
      util::ShPtr< command::FlagStatic> m_PymolColorNodesByDescription;
      //! the filename where the descriptions for the members being clustered can be found
      util::ShPtr< command::ParameterInterface> m_DescriptionsFilename;
      //! the minimium value that should be seen out of all the descriptions of the members being clustered
      util::ShPtr< command::ParameterInterface> m_MinimumDescriptor;
      //! the maximium value that should be seen out of all the descriptions of the members being clustered
      util::ShPtr< command::ParameterInterface> m_MaximumDescriptor;

      //! specifies that Nodes with less than a certain number of members will not be displayed
      util::ShPtr< command::FlagStatic> m_PymolRemoveNodesBelowSize;
      //! the minimum number of members a node must have in order to be displayed
      util::ShPtr< command::ParameterInterface> m_MinimumNumberMembers;

      //! specifies that nodes which are internally similar will not be displayed
      //! ( e.g. nodes with RMSD less than 2.3; uses binary predicate to determine comparison)
      util::ShPtr< command::FlagStatic> m_PymolRemoveInternallySimilarNodes;
      //! the girth which is used for the cutoff of similarity
      //! the binary predicate determines if nodes above or below are displayed
      util::ShPtr< command::ParameterInterface> m_SimilarityCutoff;

      //! specifies that nodes which are internally dissimilar (unclustered) will not be displayed
      //! ( e.g. nodes with RMSD above than 5.0; uses binary predicate to determine comparison)
      util::ShPtr< command::FlagStatic> m_PymolRemoveInternallyDissimilarNodes;
      //! the girth which is used for the cutoff of dissimilarity
      //! the binary predicate determines if nodes above or below are displayed
      util::ShPtr< command::ParameterInterface> m_DissimilarityCutoff;

      //! whether to stagger nodes that have no peer branch, e.g. if the other member was pruned
      util::ShPtr< command::FlagStatic> m_StaggerNodesWithMissingBranchFlag;

      //! specifies that for visualization, the clustered strings are pdbs and can be displayed as protein models
      util::ShPtr< command::FlagStatic> m_LabelProteinModelFromString;
      //! the prefix that must be added to the members in order to access the pdb file
      util::ShPtr< command::ParameterInterface> m_ProteinModelFromStringPrefix;
      //! the postfix that must be added to the members in order to access the pdb file
      util::ShPtr< command::ParameterInterface> m_ProteinModelFromStringPostfix;

      //! specifies that for visualization, the clustered strings are molecules and can be displayed from sdf file
      util::ShPtr< command::FlagStatic> m_LabelSmallMolecule;
      //! the sdf file where the molecules are going to come from
      util::ShPtr< command::ParameterInterface> m_SdfInputFile;

      //! specifies for visualization that the strings used for clustering should be used to label nodes
      util::ShPtr< command::FlagStatic> m_LabelString;

      //! specifies for visualization that the radius of the cylinder depicting a node should scale with the number
      //! of members in that node
      util::ShPtr< command::FlagInterface> m_PymolScaleNodeWithSize;

      //! bool indicating if all the members of a labeled node should be listed in a file
      util::ShPtr< command::FlagStatic> m_FlagOutputNodeMembers;
      util::ShPtr< command::ParameterInterface> m_OutputNodeMembersOutFile;

      // flag to set the minimum and maximum girth that the dendrogram in pymol reaches to - visualization only flag
      util::ShPtr< command::FlagStatic> m_SetMinMaxGirth;
      //! the minimum girth the dendrogram will go to
      util::ShPtr< command::ParameterInterface> m_MinGirth;
      //! the maximum girth the dendrogram will go to
      util::ShPtr< command::ParameterInterface> m_MaxGirth;

      //! flag indicating that preclustering should be done to speed up clustering
      util::ShPtr< command::FlagStatic> m_Precluster;

      //! the threshold that should be used to combine nodes during the preclustering step
      util::ShPtr< command::ParameterInterface> m_PreclusteringThreshold;

    public:

      // instantiate enumerator for GenerateDataset class
      static const ApplicationType Cluster_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Cluster();

      //! @brief Clone function
      //! @return pointer to new Cluster
      Cluster *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief OutputResults outputs all the results of the clustering as specified by the output formats specified
      //!        by the command line
      //! @param NODES a storage::List which contains all the nodes in the dendrogram
      //! @param OBJECT_DISTANCE_CALCULATOR the method that calculates distances between clustered objects
      //! @param COMP_FUNCTION the comparison function used to compare distances between clustered objects
      void OutputResults
      (
        const util::SiPtrList< const cluster::Node< std::string, float> > &NODES,
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, float>
        > &OBJECT_DISTANCE_CALCULATOR,
        const math::Comparisons< float>::Comparison &COMP_FUNCTION
      ) const;

      //! @brief OutputPymol outputs the dendrogram in a Python script that can be displayed by Pymol
      //! @param NODE a the Node which will be output and contains all the nodes in the dendrogram
      //! @param OBJECT_DISTANCE_CALCULATOR the method that calculates distances between clustered objects
      //! @param COMPARISON_FUNCTION the comparison function used to compare distances between clustered objects
      void OutputPymol
      (
        const cluster::Node< std::string, float> &NODE,
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, float>
        > &OBJECT_DISTANCE_CALCULATOR,
        const util::ShPtr< util::BinaryFunctionInterface< float, float, bool> > &COMPARISON_FUNCTION
      ) const;

      //! @brief GetNodeColorer provides the node colorer created according to the user specifications
      //! @return returns a ShPtr to a function interface which takes a node and returns three color points
      util::ShPtr
      <
        util::FunctionInterface< cluster::Node< std::string, float>, linal::Vector3D>
      > GetNodeColorer() const;

    }; // class Cluster

  } // namespace app
} // namespace bcl

#endif // BCL_APP_CLUSTER_H_
