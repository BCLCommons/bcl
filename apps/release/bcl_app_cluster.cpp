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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_app_cluster.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_dendrogram.h"
#include "cluster/bcl_cluster_input_classes.h"
#include "cluster/bcl_cluster_linkage_classes.h"
#include "cluster/bcl_cluster_node_colorer.h"
#include "cluster/bcl_cluster_node_description_average.h"
#include "cluster/bcl_cluster_node_description_from_file.h"
#include "cluster/bcl_cluster_output_classes.h"
#include "cluster/bcl_cluster_output_pymol.h"
#include "cluster/bcl_cluster_output_pymol_label_protein_model_from_string.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_color_gradient.h"

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

    const ApplicationType Cluster::Cluster_Instance
    (
      GetAppGroups().AddAppToGroup( new Cluster(), GetAppGroups().e_Bcl)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // default constructor
    Cluster::Cluster() :
      m_InputFile
      (
        new command::FlagStatic
        (
          "distance_input_file", "This flag is for using a file of distances as input to the clustering",
           command::Parameter
          (
            "distance_input_file",
            "The file which contains the distances between the objects to be clustered"
          )
        )
      ),
      m_InputFormat
      (
        new command::FlagStatic
        (
          "input_format",
          "This flag is for specifying the format of the distance_input_file",
            command::Parameter
            (
              "format",
              "the format of the input data",
              command::ParameterCheckEnumerate< cluster::InputClasses< std::string, float> >(),
              cluster::GetInputClasses< std::string, float>().e_TableLowerTriangle
            )
         )
      ),
      m_OutputFormat
      (
        new command::FlagDynamic
        (
          "output_format",
          "The format(s)  in which that the dendrogram/clusters should be output",
           command::Parameter
          (
            "format_list",
            "Give a list of the formats in which the dendrogram/clusters should be output",
            command::ParameterCheckEnumerate< cluster::OutputClasses< std::string, float> >(),
            cluster::GetOutputClasses< std::string, float>().e_Rows
          ), 1, cluster::GetOutputClasses< std::string, float>().GetEnumCount()
        )
      ),
      m_HeightCutoff
      (
        new command::FlagStatic
        (
          "height_cutoff", "This flag is for stopping the clustering once the has reached a certain cluster girth",
          command::Parameter
          (
            "height_cutoff_value",
            "The value which is desired for the height cutoff - default undefined i.e. no height cutoff",
            util::Format()( util::GetUndefined< float>())
          )
        )
      ),
      m_OutputFile
      (
        new command::FlagStatic
        (
          "output_file",
          "Path and name of the output file which will hold the results of the clustering. The ouput format is appended to the file name.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "cluster.txt"
          )
        )
      ),
      m_Linkage
      (
        new command::FlagStatic
        (
          "linkage", "This flag is for specifying the type of linkage to be used during clustering",
           command::Parameter
          (
            "linkage_type",
            "The type of linkage to use during clustering.",
            command::ParameterCheckEnumerate< cluster::LinkageClasses< std::string, float> >(),
            cluster::GetLinkageClasses< std::string, float>().e_Single
          )
        )
      ),
      m_DistanceType
      (
        new command::FlagStatic
        (
          "distance_definition",
          "Defines if the distance measure is a similarity or dissimilarity measure.",
          command::Parameter
          (
            "distance measure type",
            "If two objects with a large distance are more similar than two objects with a small distance, then \"greater\" should be used. If two objects with a small distance are more similar than two objects with a large distance, then \"less\" should be used.",
            command::ParameterCheckEnumerate< math::Comparisons< float> >(),
            math::Comparisons< float>::GetEnums().e_Less.GetName()
          )
        )
      ),
      m_OutputPymol
      (
        new command::FlagStatic
        (
          "output_pymol",
          "output the dendrogram in a Python script that can be displayed by Pymol"
        )
      ),
      m_CylinderLength
      (
        new command::Parameter
        (
          "cylinder_length",
          "the unit length of a cylinder used to build the dendrogram",
          "100"
        )
      ),
      m_CylinderRadius
      (
        new command::Parameter
        (
          "cylinder_radius",
          "the unit radius of a cylinder used to build the dendrogram",
          "25"
        )
      ),
      m_CylinderSeparation
      (
        new command::Parameter
        (
          "cylinder_separation",
          "the unit separation between cylinders used to build the dendrogram",
          "50"
        )
      ),
      m_MaxNodeLabels
      (
        new command::Parameter
        (
          "max_number_node_labels",
          "the maximum number of labels you want to be displayed in your dendrogram",
          "100"
        )
      ),
      m_NumberLabelsAlongYScale
      (
        new command::Parameter
        (
          "number_labels_along_y_scale",
          "number of labels you want to be placed along the vertical axis of your dendrogram indicating node heigh levels",
          "10"
        )
      ),
      m_OutputPymolOutputFile
      (
        new command::Parameter
        (
          "pymol_ouput_file",
          "the name of the file which will have the Python code for displaying the dendrogram in Pymol",
          "dendrogram.py"
        )
      ),
      m_PymolColorNodesByDescription
      (
        new command::FlagStatic
        (
          "pymol_color_nodes_by_description",
          "for visualization, nodes will be colored according to some numerical description provided by a file"
        )
      ),
      m_DescriptionsFilename
      (
        new command::Parameter
        (
          "descriptions_filename",
          "the file which contains the descriptions for the objects being clustered. Descriptions must be positive. Color scheme goes red(small)->yellow(middle)->white(large) description value. Range from small to large is given by parameters below. Format per line: <member> <description>",
          "descriptions.txt"
        )
      ),
      m_MinimumDescriptor
      (
        new command::Parameter
        (
          "minimum_descriptor",
          "the minimium value that should be seen out of all the descriptions of the members being clustered",
          "0.0"
        )
      ),
      m_MaximumDescriptor
      (
        new command::Parameter
        (
          "maximum_descriptor",
          "the maximium value that should be seen out of all the descriptions of the members being clustered",
          "100000.0"
        )
      ),
      m_PymolRemoveNodesBelowSize
      (
        new command::FlagStatic
        (
          "remove_nodes_below_size",
          "nodes with less than a certain number of members will be removed"
        )
      ),
      m_MinimumNumberMembers
      (
        new command::Parameter
        (
          "minimum_number_members",
          "The minimum number of members a node must have to be kept",
          "1"
        )
      ),
      m_PymolRemoveInternallySimilarNodes
      (
        new command::FlagStatic
        (
          "remove_internally_similar_nodes",
          "nodes whose members are too similar will be removed"
        )
      ),
      m_SimilarityCutoff
      (
        new command::Parameter
        (
          "similarity_cutoff",
          "the girth which is used for the cutoff of similarity; whether nodes above or below the cutoff are removed is automatically handled by the binary predicate specified with the binary predicate option above",
          util::Format()( util::GetUndefined< float>())
        )
      ),
      m_PymolRemoveInternallyDissimilarNodes
      (
        new command::FlagStatic
        (
          "remove_internally_dissimilar_nodes",
          "allows removing nodes that do not form a logical cluster due to poor linkage "
          "When pruning large dendrograms with all 3 removal flags, this flag often works well with with a value just "
          "below the value given to -remove_internally_similar_nodes (if distance_type is a dissimilarity measure, slightly above)"
        )
      ),
      m_DissimilarityCutoff
      (
        new command::Parameter
        (
          "dissimilarity_cutoff",
          "the girth which is used for the cutoff of dissimilarity; whether nodes above or below the cutoff are removed is automatically handled by the binary predicate specified with the binary predicate option above",
          util::Format()( util::GetUndefined< float>())
        )
      ),
      m_StaggerNodesWithMissingBranchFlag
      (
        new command::FlagStatic
        (
          "stagger_pruned_branches",
          "whether to stagger nodes that have no peer branch, e.g. if the other member was pruned"
        )
      ),
      m_LabelProteinModelFromString
      (
        new command::FlagStatic
        (
          "pymol_label_output_protein_model_from_string",
          "output the protein models from the strings in the table"
        )
      ),
      m_ProteinModelFromStringPrefix
      (
        new command::Parameter
        (
          "model_from_string_prefix",
          "the prefix that must be added to the members in order to access the file",
          ""
        )
      ),
      m_ProteinModelFromStringPostfix
      (
        new command::Parameter
        (
          "model_from_string_postfix",
          "the postfix that must be added to the members in order to access the file",
          ""
        )
      ),
      m_LabelSmallMolecule
      (
        new command::FlagStatic
        (
          "pymol_label_output_small_molecule",
          "label nodes with small molecules. The rows and columns of the input distance file is assumed to have the molecules numbered from "
          "0 to N-1, where N is the number of molecules and this numbering should correspond to the order in which the molecules "
          "are found in the sdf file."
        )
      ),
      m_SdfInputFile
      (
        new command::Parameter
        (
          "sdf_input_file",
          "the sdf file from which the molecules will be taken",
          "molecules.sdf"
        )
      ),
      m_LabelString
      (
        new command::FlagStatic
        (
          "pymol_label_output_string",
          "output labels which are just the strings in the table header. This is the default labeler if no other label_output_pymol flags are given."
        )
      ),
      m_PymolScaleNodeWithSize
      (
        new command::FlagStatic
        (
          "pymol_scale_node_with_size",
          "scale the diameter of the node with the number of members which are contained in that node"
        )
      ),
      m_FlagOutputNodeMembers
      (
        new command::FlagStatic
        (
          "output_node_members",
          "if set then labeled nodes will have members printed out list files with filename node_<node_identifier>.ls. This will automatically create pymol output. The number of files output can be adjusted by changing the max number of node labels parameter in the output_pymol flag."
        )
      ),
      m_OutputNodeMembersOutFile
      (
        new command::Parameter
        (
          "output_node_members_outfile",
          "the base filename for the files that will contain the members of each node. If pymol output is also specified, the base filename is the same as the user provided dendrogram filename.",
          "output_members_file.txt"
        )
      ),
      m_SetMinMaxGirth
      (
        new command::FlagStatic
        (
          "pymol_set_min_max_girth",
          "define the minimum and maximum girth that is labeled along the y-axis. The dendrogram will span this range"
        )
      ),
      m_MinGirth
      (
        new command::Parameter
        (
          "minimum_girth",
          "the minimum girth value that will be labeled on the y-axis",
          util::Format()( util::GetUndefined< double>())
        )
      ),
      m_MaxGirth
      (
        new command::Parameter
        (
          "maximum_girth",
          "the maximum girth value that will be labeled on the y-axis",
          util::Format()( util::GetUndefined< double>())
        )
      ),
      m_Precluster
      (
        new command::FlagStatic
        (
          "precluster",
          "before clustering is done, do a preclustering step to reduce the number of iterations that must occur"
        )
      ),
      m_PreclusteringThreshold
      (
        new command::Parameter
        (
          "threshold",
          "how similar objects must be in order to be preclustered together",
          util::Format()( util::GetUndefined< float>())
        )
      )
    {
      // add parameter to "m_PymolRemoveNodesBelowSize"
      m_PymolRemoveNodesBelowSize->PushBack( m_MinimumNumberMembers);

      // add parameter to "m_PymolRemoveInternallySimilarNodes"
      m_PymolRemoveInternallySimilarNodes->PushBack( m_SimilarityCutoff);
      m_PymolRemoveInternallyDissimilarNodes->PushBack( m_DissimilarityCutoff);

      // add the parameters to "m_OutputPymol"
      m_OutputPymol->PushBack( m_CylinderLength);
      m_OutputPymol->PushBack( m_CylinderRadius);
      m_OutputPymol->PushBack( m_CylinderSeparation);
      m_OutputPymol->PushBack( m_MaxNodeLabels);
      m_OutputPymol->PushBack( m_NumberLabelsAlongYScale);
      m_OutputPymol->PushBack( m_OutputPymolOutputFile);

      // add parameters to "m_PymolColorNodesByDescription"
      m_PymolColorNodesByDescription->PushBack( m_DescriptionsFilename);
      m_PymolColorNodesByDescription->PushBack( m_MinimumDescriptor);
      m_PymolColorNodesByDescription->PushBack( m_MaximumDescriptor);

      // add parameters to "m_LabelProteinModelFromString"
      m_LabelProteinModelFromString->PushBack( m_ProteinModelFromStringPrefix);
      m_LabelProteinModelFromString->PushBack( m_ProteinModelFromStringPostfix);

      // add parameters
      m_LabelSmallMolecule->PushBack( m_SdfInputFile);

      // add flag
      m_FlagOutputNodeMembers->PushBack( m_OutputNodeMembersOutFile);

      // add parameters to "m_SetMinMaxGirth"
      m_SetMinMaxGirth->PushBack( m_MinGirth);
      m_SetMinMaxGirth->PushBack( m_MaxGirth);

      // the threshold that should be used to combine nodes during the preclustering step
      m_Precluster->PushBack( m_PreclusteringThreshold);
    }

    //! @brief Clone function
    //! @return pointer to new Quality
    Cluster *Cluster::Clone() const
    {
      return new Cluster( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Cluster::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string Cluster::GetDescription() const
    {
      return "Hierarchical agglomerative clustering of data";
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> Cluster::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input file flag
      sp_cmd->AddFlag( m_InputFile);

      // input file format flag
      sp_cmd->AddFlag( m_InputFormat);

      // output format flag
      sp_cmd->AddFlag( m_OutputFormat);

      // input file format flag
      sp_cmd->AddFlag( m_Linkage);

      // height cutoff flag
      sp_cmd->AddFlag( m_HeightCutoff);

      // height output file flag
      sp_cmd->AddFlag( m_OutputFile);

      // distance type flag
      sp_cmd->AddFlag( m_DistanceType);

      // add "m_PymolRemoveNodesBelowSize" to command object
      sp_cmd->AddFlag( m_PymolRemoveNodesBelowSize);

      // add "m_PymolRemoveInternallySimilarNodes" to command object
      sp_cmd->AddFlag( m_PymolRemoveInternallySimilarNodes);
      sp_cmd->AddFlag( m_PymolRemoveInternallyDissimilarNodes);

      sp_cmd->AddFlag( m_StaggerNodesWithMissingBranchFlag);

      // add "m_OutputPymol" to the command
      sp_cmd->AddFlag( m_OutputPymol);

      // add "m_PymolColorNodesByDescription" to command object
      sp_cmd->AddFlag( m_PymolColorNodesByDescription);

      sp_cmd->AddFlag( m_LabelProteinModelFromString);

      sp_cmd->AddFlag( m_LabelSmallMolecule);

      // add "m_LabelProteinModelFromString" to command object
      sp_cmd->AddFlag( m_LabelString);

      // add "m_PymolScaleNodeWithSize" to command object
      sp_cmd->AddFlag( m_PymolScaleNodeWithSize);

      sp_cmd->AddFlag( m_FlagOutputNodeMembers);

      // add "m_SetMinMaxGirth" to command object
      sp_cmd->AddFlag( m_SetMinMaxGirth);

      // flag indicating that preclustering should be done to speed up clustering
      sp_cmd->AddFlag( m_Precluster);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &Cluster::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::Cluster, terms of use, appropriate citation, installation "
        "procedures, BCL::Cluster execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::Cluster?\n"
        "BCL::Cluster is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::Cluster is a clustering utility. It reads in "
        "pairwise distances for a set of objects from a file.  It then uses a hierarchical agglomerative clustering "
        "algorithm to cluster the dataset.  Pairwise distances can be similarity or dissimilarity measures.\n"
        "BCL::Cluster outputs information about the clustering hierarchy in text files.  In addition, BCL::Cluster "
        "outputs a file which can be used in conjunction with the Pymol Molecular Graphics System in order to "
        "visualize the resulting dendrogram.  As an option, the actual biological molecules that were being clustered "
        "(if biological molecules were being clustered) can be visualized with the dendrogram.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Cluster.\n"
        "When using BCL::Cluster in a publication, please cite the following publication describing the application's "
        "development:\n"
        "\n"
        "N. Alexander, N. Woetzel, and J. Meiler, Bcl::Cluster: A method for clustering biological molecules coupled "
        "with visualization in the Pymol Molecular Graphics System, in 2011 IEEE 1st International Conference on "
        "Computational Advances in Bio and Medical Sciences (ICCABS), 2011, pp. 13–18."
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::Cluster.\n"
        "Running BCL::Cluster consists of two main steps.\n"
        "\n"
        "1) Create a file containing the pairwise distances between objects to be clustered.\n"
        "\n"
        "2) Run BCL::Cluster to cluster on the pairwise distance file\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::Cluster executable program.\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl_cluster.exe -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        "For more general information about the product, type bcl_cluster.exe -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::Cluster.\n"
        "BCL::Cluster is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE"
        "<bcl.exe>  Cluster -distance_input_file models_distance.matrix -input_format TableLowerTriangle -output_file "
        "models.out -output_format Rows"
        + DefaultSectionSeparator()
      );
      return readme;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int Cluster::Main() const
    {
      // create const string "input_format" and initialize with the input format from the command line
      const std::string input_format( m_InputFormat->GetFirstParameter()->GetValue());

      BCL_MessageDbg( "input format is " + input_format);

      // create ShPtr to an cluster::InputInterface "input_interface"
      util::ShPtr< cluster::InputInterface< std::string, float> > input_interface
      (
        // initialize with the InputClass given by "input_format"
        ( *cluster::GetInputClasses< std::string, float>().GetEnumFromName( input_format))->Clone()
      );

      // make sure that "input_interface" is defined i.e. an input class was succesfully gotten
      BCL_Assert( input_interface.IsDefined(), "the InputInterface ShPtr is not defined");

      // create string "input_filename" and initialize with the name of the input file from the command line
      const std::string input_filename( m_InputFile->GetFirstParameter()->GetValue());
      // create ShPtr to a linkage interface "linkage" and initialize with the Linkage given by the command line

      util::ShPtr< cluster::LinkageInterface< std::string, float> > linkage
      (
        (
          *cluster::GetLinkageClasses< std::string, float>().GetEnumFromName
          (
            m_Linkage->GetFirstParameter()->GetValue()
          )
        )->Clone()
      );

      // make sure that "linkage" is defined
      BCL_Assert( linkage.IsDefined(), "the linkage ShPtr is not defined");

      BCL_MessageCrt( "start getting distance function");
      // create IFStream "ifstream"
      io::IFStream ifstream;
      io::File::MustOpenIFStream( ifstream, input_filename);

      // create ShPtr to a function interface "object_distance_calculator" and initialize with a DistancesStored object
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, float>
      > object_distance_calculator
      (
        // initialize with the distances data returned by "input_interface"
        input_interface->HandleInput( ifstream)
      );
      io::File::CloseClearFStream( ifstream);

      // make sure that "object_distance_calculator" is defined
      BCL_Assert( object_distance_calculator.IsDefined(), "the object_distance_calculator ShPtr is not defined");

      // set the DistanceFunction data member of "linkage", which calculates the distance between any two
      // objects being clustered
      linkage->SetDistanceFunction( object_distance_calculator);

      // create const float "height_cutoff" and initialize with the height cutoff given by the command line
      float height_cutoff( m_HeightCutoff->GetFirstParameter()->GetNumericalValue< float>());

      // create math::Comparisons "comp_function" and initialize with command line in order to define
      // the comparison function i.e. are similarities or dissimilarities used in the clustering
      // "<" is used for dissimilarities ( e.x. RMSD) and ">" is used for similarities
      const math::Comparisons< float>::Comparison comp_function( m_DistanceType->GetFirstParameter()->GetValue());
      const util::BinaryFunctionInterface< float, float, bool> &comp_girth_fun( **comp_function);

      BCL_MessageDbg( "the comparison operator will be " + util::Format()( comp_function) + " it is defined " + util::Format()( util::IsDefined( height_cutoff)));

      // true if the height cutoff has not been set by user and the comparison function has been set to indicate the
      // distance measures are similarity measures.
      // need to set the height cutoff to a very negative number so that clustering will not stop due to height cutoff
      if
      (
        !util::IsDefined( height_cutoff) &&
        (
          comp_function == math::Comparisons< float>::GetEnums().e_Greater ||
          comp_function == math::Comparisons< float>::GetEnums().e_GreaterEqual
        )
      )
      {
        BCL_MessageStd
        (
          "height cutoff undefined with greater or greater_equal comparison so setting to most negative float possible"
        );
        // set height cutoff to the most negative number possible for a float
        height_cutoff = -std::numeric_limits< float>::max();
      }
      else if
      (
        !util::IsDefined( height_cutoff) &&
        (
          comp_function == math::Comparisons< float>::GetEnums().e_Less ||
          comp_function == math::Comparisons< float>::GetEnums().e_LessEqual
        )
      )
      {
        BCL_MessageStd
        (
          "height cutoff undefined with less or less_equal comparison so setting to most positive float possible"
        );
        // set height cutoff to the most positive number possible for a float
        height_cutoff = std::numeric_limits< float>::max();
      }

      // set the comparison operator in the linkage calculator
      linkage->SetBinaryPredicate( comp_girth_fun);

      // get the threshold for combining nodes during preclustering
      const float preclustering_threshold( m_PreclusteringThreshold->GetNumericalValue< float>());

      // create util::Dendrogram "dendrogram"; after the construction, the clustering will be complete
      const cluster::Dendrogram< std::string, float> dendrogram
      (
        linkage, //< the linkage to be used
        input_interface->GetInputObjects(), //< the actual objects to be clustered
        height_cutoff, //< height to stop clustering
        **comp_function, //< comparison function to be used
        preclustering_threshold, //< threshold for preclustering
        m_SimilarityCutoff->GetNumericalValue< float>() //< threshold for removing internally similar nodes
      );

      // output the girth of the top node
      BCL_MessageDbg( "highest girth is " + util::Format()( dendrogram.GetNode().GetGirth()));

      cluster::Node< std::string, float> node( dendrogram.GetNode());

      // true if it is desired to not show nodes which have a small number of members
      if( m_PymolRemoveNodesBelowSize->GetFlag())
      {
        node.RemoveNodesBelowSize( m_MinimumNumberMembers->GetNumericalValue< size_t>());
      }

      // true if it is desired to not show leaf nodes that do not represent a sufficiently well-linked cluster
      if( m_PymolRemoveInternallyDissimilarNodes->GetFlag())
      {
        if( !m_PymolRemoveNodesBelowSize->GetFlag() && !m_PymolRemoveInternallySimilarNodes->GetFlag())
        {
          BCL_MessageCrt
          (
            m_PymolRemoveInternallyDissimilarNodes->GetName() + " has no effect unless "
            + m_PymolRemoveNodesBelowSize->GetName() + " or "
            + m_PymolRemoveInternallySimilarNodes->GetName() + " is used"
          );
        }
        else
        {
          BCL_Assert
          (
            !m_StaggerNodesWithMissingBranchFlag->GetFlag(),
            m_PymolRemoveInternallyDissimilarNodes->GetName() + " is incompatible with "
            + m_StaggerNodesWithMissingBranchFlag->GetName()
          );
          node.RemoveSingularBranches();
          node.RemoveNodesWithLowSimilarity( m_DissimilarityCutoff->GetNumericalValue< float>(), comp_girth_fun);
        }
      }
      if( !m_StaggerNodesWithMissingBranchFlag->GetFlag())
      {
        node.RemoveSingularBranches();
      }

      const util::SiPtrList< const cluster::Node< std::string, float> > node_list
      (
        1, util::SiPtr< const cluster::Node< std::string, float> >( node)
      );

      // Output the results of the clustering
      OutputResults( node_list, object_distance_calculator, comp_function);

      BCL_MessageDbg( "finished outputting to formats ");

      if( m_OutputPymol->GetFlag() || m_FlagOutputNodeMembers->GetFlag())
      {
        OutputPymol( node, object_distance_calculator, *comp_function);
      }

      //successful end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Cluster::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Cluster::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief OutputResults outputs all the results of the clustering as specified by the output formats specified
    //!        by the command line
    //! @param NODES a storage::List which contains all the nodes in the dendrogram
    //! @param OBJECT_DISTANCE_CALCULATOR the method that calculates distances between clustered objects
    //! @param COMP_FUNCTION the comparison function used to compare distances between clustered objects
    void Cluster::OutputResults
    (
      const util::SiPtrList< const cluster::Node< std::string, float> > &NODES,
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, float>
      > &OBJECT_DISTANCE_CALCULATOR,
      const math::Comparisons< float>::Comparison &COMP_FUNCTION
    ) const
    {
      // create string "output_file_name" and initialize with the command line argument "m_OutputFile"
      const std::string output_file_name( m_OutputFile->GetFirstParameter()->GetValue());

      // create string "output_file_name_without_extension" and initialize with an extensionless "output_file_name"
      const std::string output_file_name_without_extension( io::File::RemoveFullExtension( output_file_name));

      // create string "output_file_name_extension" and initialize with an extension of "output_file_name"
      const std::string output_file_name_extension( io::File::GetFullExtension( output_file_name));

      // loop over the output formats given in the command line so that each format can be outputted
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          itr( m_OutputFormat->GetParameterList().Begin()), itr_end( m_OutputFormat->GetParameterList().End());
        itr != itr_end;
        ++itr
      )
      {
        // create an OutputClass "current_output_format"
        cluster::OutputClasses< std::string, float>::OutputClass current_output_format
        (
          // initialize with the desired output format indicated by "itr"
          cluster::GetOutputClasses< std::string, float>().GetEnumFromName( ( *itr)->GetValue())
        );

        if( current_output_format == cluster::GetOutputClasses< std::string, float>().e_Centers)
        {
          cluster::OutputCenters< std::string, float> output_center
          (
            OBJECT_DISTANCE_CALCULATOR, *COMP_FUNCTION
          );

          // create const string "temp_output_file_name" and initialize with "output_file_name_without_extension"
          // the name of "current_output_format" and "output_file_name_extension"
          const std::string temp_output_file_name
          (
            output_file_name_without_extension + ".Centers." + output_file_name_extension
          );

          output_center.WriteOutput( temp_output_file_name, NODES);

          continue;
        }
        else if
        (
          current_output_format == cluster::GetOutputClasses< std::string, float>().e_Matrix
          || current_output_format == cluster::GetOutputClasses< std::string, float>().e_Table
        )
        {
          cluster::OutputSortedMatrix< std::string, float> output_center
          (
            OBJECT_DISTANCE_CALCULATOR,
            *COMP_FUNCTION,
            current_output_format == cluster::GetOutputClasses< std::string, float>().e_Table
          );

          // create const string "temp_output_file_name" and initialize with "output_file_name_without_extension"
          // the name of "current_output_format" and "output_file_name_extension"
          const std::string temp_output_file_name
          (
            output_file_name_without_extension + "." + current_output_format.GetName() + "." + output_file_name_extension
          );

          output_center.WriteOutput( temp_output_file_name, NODES);

          continue;
        }

        // create const string "temp_output_file_name" and initialize with "output_file_name_without_extension"
        // the name of "current_output_format" and "output_file_name_extension"
        const std::string temp_output_file_name
        (
          output_file_name_without_extension + "." + current_output_format.GetName() + "." + output_file_name_extension
        );

        // output the filename the current format is being written to
        BCL_MessageStd( "Outputting to file " + temp_output_file_name);

        // create OFStream "write"
        io::OFStream write;

        // write to "temp_output_file_name"
        ( *current_output_format)->WriteOutput( temp_output_file_name, NODES);
      }
    }

    //! @brief OutputPymol outputs the dendrogram in a Python script that can be displayed by Pymol
    //! @param NODE a the Node which will be output and contains all the nodes in the dendrogram
    //! @param OBJECT_DISTANCE_CALCULATOR the method that calculates distances between clustered objects
    //! @param COMPARISON_FUNCTION the comparison function used to compare distances between clustered objects
    void Cluster::OutputPymol
    (
      const cluster::Node< std::string, float> &NODE,
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, float>
      > &OBJECT_DISTANCE_CALCULATOR,
      const util::ShPtr< util::BinaryFunctionInterface< float, float, bool> > &COMPARISON_FUNCTION
    ) const
    {
      // create method for coloring node NodeColorer "node_colorer" and initialize with GetNodeColorer function
      util::ShPtr
      <
        util::FunctionInterface< cluster::Node< std::string, float>, linal::Vector3D>
      > node_colorer( GetNodeColorer());

      // create method for creating a label for a node OutputPymolLabelProteinModelFromString "label_maker"
      // and initialize with the OutputPymolLabelString as a default, but could be changed below
      // Right now this gives the same behavior if m_LabelString is passed, but by putting a default value here
      // it ensures that the label maker is defined incase the user does not activate one of the label maker flags
      util::ShPtr
      <
        util::FunctionInterface
        <
          storage::Triplet //!< argument type
          <
            util::SiPtr< const cluster:: Node< std::string, float> >,
            storage::VectorND< 4, double>,
            std::string
          >,
          std::string //!< return type
        >
      > label_maker
      (
        util::ShPtr
        <
          util::FunctionInterface
          <
            storage::Triplet //!< argument type
            <
              util::SiPtr< const cluster:: Node< std::string, float> >,
              storage::VectorND< 4, double>,
              std::string
            >,
            std::string //!< return type
          >
        >( new cluster::OutputPymolLabelString< float>( OBJECT_DISTANCE_CALCULATOR, COMPARISON_FUNCTION, util::GetColors().e_Magenta, m_FlagOutputNodeMembers->GetFlag()))
      );

      // if the label maker should actually be a ProteinModelFromString label maker then set "label_maker" accordingly
      if( m_LabelProteinModelFromString->GetFlag())
      {
        label_maker =
        util::ShPtr
        <
          util::FunctionInterface
          <
            storage::Triplet //!< argument type
            <
              util::SiPtr< const cluster:: Node< std::string, float> >,
              storage::VectorND< 4, double>,
              std::string
            >,
            std::string //!< return type
          >
        >
        (
          new cluster::OutputPymolLabelProteinModelFromString< float>
          (
            OBJECT_DISTANCE_CALCULATOR,     COMPARISON_FUNCTION,
            m_ProteinModelFromStringPrefix->GetValue(), m_ProteinModelFromStringPostfix->GetValue(),
            util::GetColors().e_Magenta, m_FlagOutputNodeMembers->GetFlag()
          )
        );
      }
      // true if the label maker should be a small molecule type
      else if( m_LabelSmallMolecule->GetFlag())
      {
        label_maker =
        util::ShPtr
        <
          util::FunctionInterface
          <
            storage::Triplet //!< argument type
            <
              util::SiPtr< const cluster:: Node< std::string, float> >,
              storage::VectorND< 4, double>,
              std::string
            >,
            std::string //!< return type
          >
        >
        (
          new cluster::OutputPymolLabelSmallMolecule< float>
          (
            OBJECT_DISTANCE_CALCULATOR,     COMPARISON_FUNCTION,
            m_SdfInputFile->GetValue(), util::GetColors().e_Magenta,
            m_FlagOutputNodeMembers->GetFlag()
          )
        );
      }
      // else if the label should just be a string, then set "label_maker" accordingly
      else if( m_LabelString->GetFlag())
      {
        label_maker = util::ShPtr
        <
          util::FunctionInterface
          <
            storage::Triplet //!< argument type
            <
              util::SiPtr< const cluster:: Node< std::string, float> >,
              storage::VectorND< 4, double>,
              std::string
            >,
            std::string //!< return type
          >
        >( new cluster::OutputPymolLabelString< float>( OBJECT_DISTANCE_CALCULATOR, COMPARISON_FUNCTION, util::GetColors().e_Magenta, m_FlagOutputNodeMembers->GetFlag()));
      }

      bool scale_with_node_size( m_PymolScaleNodeWithSize->GetFlag());

      const storage::VectorND< 2, float> min_max_girth
      (
        m_MinGirth->GetNumericalValue< float>(), m_MaxGirth->GetNumericalValue< float>()
      );

      // create method for writing the Python script necessary for displaying the dendrogram in Pymol
      // initialize with the parameters from the command line and "node_colorer" and "label_maker"
      cluster::OutputPymol< std::string, float> output_pymol
      (
        m_CylinderLength->GetNumericalValue< double>(),         m_CylinderRadius->GetNumericalValue< double>(),
        m_CylinderSeparation->GetNumericalValue< double>(),     node_colorer,
        m_MaxNodeLabels->GetNumericalValue< size_t>(),          label_maker,
        util::GetColors().e_Magenta,                            m_NumberLabelsAlongYScale->GetNumericalValue< size_t>(),
        COMPARISON_FUNCTION->operator()( 0, 1),                 min_max_girth, scale_with_node_size
      );

      // write the Python script necessary to display "NODE" in Pymol into "m_OutputPymolOutputFile"
      std::string pymol_output_filename( m_OutputPymolOutputFile->GetValue());
      if( m_FlagOutputNodeMembers->GetFlag() && !m_OutputPymol->GetFlag())
      {
        const std::string output_member_filename( m_OutputNodeMembersOutFile->GetValue());
        pymol_output_filename = output_member_filename;
      }

      output_pymol.WriteOutput
      (
        pymol_output_filename,
        util::SiPtrList< const cluster::Node< std::string, float> >
        (
          1, util::SiPtr< const cluster::Node< std::string, float> >( NODE)
        )
      );
    }

    //! @brief GetNodeColorer provides the node colorer created according to the user specifications
    //! @return returns a ShPtr to a function interface which takes a node and returns three color points
    util::ShPtr
    <
      util::FunctionInterface< cluster::Node< std::string, float>, linal::Vector3D>
    > Cluster::GetNodeColorer() const
    {
      // true if user desires for the nodes to be colored based on descriptions in a file provided
      if( m_PymolColorNodesByDescription->GetFlag())
      {
        // create the vector of colors to use for gradient
        const storage::Vector< util::Color> gradient_points
        (
          storage::Vector< util::Color>::Create
          (
            util::GetColors().e_Red, util::GetColors().e_Yellow, util::GetColors().e_White
          )
        );

        // create coloring function "color_function" initialize and initialize with a OutputPymolColorGradient
        const util::ShPtr< util::FunctionInterface< double, linal::Vector3D> > color_function
        (
          new util::ColorGradient
          (
            math::Range< double>
            (
              m_MinimumDescriptor->GetNumericalValue< double>(),
              m_MaximumDescriptor->GetNumericalValue< double>()
            ),
            gradient_points
          )
        );

        // create const string "descriptions_filename" and initialize with the file that has member descriptions
        const std::string descriptions_filename( m_DescriptionsFilename->GetValue());

        // create method for getting individual member descriptors
        util::ShPtr< util::FunctionInterface< std::string, double> > member_descriptor_function
        (
          new cluster::NodeDescriptionFromFile( descriptions_filename)
        );

        // create descriptor function "node_descriptor_function"
        util::ShPtr< util::FunctionInterface< cluster::Node< std::string, float>, double> >
        node_descriptor_function
        (
          new cluster::NodeDescriptionAverage< std::string, float>( member_descriptor_function)
        );

        // return method for coloring node NodeColorer initialized with
        // "node_descriptor_function" and "color_function"
        return util::ShPtr
        <
          util::FunctionInterface< cluster::Node< std::string, float>, linal::Vector3D>
        >
        (
          new cluster::NodeColorer< std::string, float>
          (
            node_descriptor_function, color_function
          )
        );
      }
      // "m_PymolColorNodesByDescription" flag is not set so just use the default colors provided by
      // "NodeColorer" in order to color the nodes
      else
      {
        // return NodeColorer method for coloring node
        return util::ShPtr
        <
          util::FunctionInterface< cluster::Node< std::string, float>, linal::Vector3D>
        >
        (
          new cluster::NodeColorer< std::string, float>()
        );
      }
    }

  } // namespace app
} // namespace bcl
