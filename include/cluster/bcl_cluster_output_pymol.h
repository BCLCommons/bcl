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

#ifndef BCL_CLUSTER_OUTPUT_PYMOL_H_
#define BCL_CLUSTER_OUTPUT_PYMOL_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_interface.h"
#include "bcl_cluster_output_pymol_label_string.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputPymol
    //! @brief is a type of Output interface which outputs the results of a Dendrogram into a Python script
    //!        The Python script can then be run by Pymol in order to display the Dendrogram.
    //!
    //! @see @link example_cluster_output_pymol.cpp @endlink
    //! @author alexanns
    //! @date August 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class OutputPymol :
      public OutputInterface< t_DataType, t_PrecisionType>
    {
    protected:

    //////////
    // data //
    //////////

      //! the unit length of the cylinders used to build the dendrogram
      double m_CylinderLength;

      //! the unit radius of the cylinders used to build the dendrogram
      double m_CylinderRadius;

      //! the diameter of the cylinders used to build the dendrogram
      double m_CylinderDiameter;

      //! the unit separation between cylinders i.e. the base leaves will be separated by this amount
      double m_CylinderSeparation;

      //! the method that is used to determine the color that a node should have
      //! the color is given by three color points (e.g. red, green, blue)
      util::ShPtr< util::FunctionInterface< Node< t_DataType, t_PrecisionType>, linal::Vector3D> > m_NodeColorer;

      //! the maximum number of labels that should be displayed on the dendrogram
      size_t m_MaxNodeLabels;

      //! the method for creating a label for nodes
      //! The Triplet argument to the util::FunctionIterface consists of the following :
      //! * util::SiPtr< const Node< t_DataType, t_PrecisionType> > is the node which will be output
      //! * storage::VectorND< 4, double> contains in this order
      //!     - the x coordinate at which the node is centered and the label will be centered
      //!     - the bottom starting y-coordinate of the node
      //!     - the top ending y-coordinate of the node
      //!     - the cylinder diameter of the node
      //! * std::string is the name of the file and path to which the the node label and all the
      //!    Python code is being written
      //! the string return type contains all the Python code necessary to display the label in Pymol
      util::ShPtr
      <
        util::FunctionInterface
        <
          storage::Triplet //!< argument type
          <
            util::SiPtr< const Node< t_DataType, t_PrecisionType> >,
            storage::VectorND< 4, double>,
            std::string
          >,
          std::string //!< return type
        >
      > m_LabelMaker;

      //! label color
      util::Color m_Color;

      //! the number of girth labels that should be displayed along the vertical axis
      size_t m_NumberLabelsAlongYScale;

      bool m_ScaleWithNodeSize;

      bool m_SmallNodesOnBottom;

      //! the minimum and maximum girth which will be displayed on the y-axis
      //! the base node will be drawn to the min or max depending on the value of "m_SmallNodesOnBottom"
      storage::VectorND< 2, t_PrecisionType> m_MinMaxGirth;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! girth label color
      static const util::Color &GetGirthLabelColor()
      {
        return util::GetColors().e_Black;
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      OutputPymol() :
        m_CylinderLength(),
        m_CylinderRadius(),
        m_CylinderDiameter(),
        m_CylinderSeparation(),
        m_NodeColorer(),
        m_MaxNodeLabels(),
        m_LabelMaker(),
        m_Color(),
        m_NumberLabelsAlongYScale(),
        m_ScaleWithNodeSize(),
        m_SmallNodesOnBottom(),
        m_MinMaxGirth( storage::VectorND< 2, t_PrecisionType>( util::GetUndefined< t_PrecisionType>(), util::GetUndefined< t_PrecisionType>()))
      {
      }

      //! @brief constructor taking all member variables
      //! @param CYLINDER_LENGTH the unit length of the cylinders used to build the dendrogram
      //! @param CYLINDER_RADIUS the radius of the cylinders used to build the dendrogram
      //! @param CYLINDER_SEPARATION the unit separation between cylinders i.e. amount separating base leaves
      //! @param NODE_COLORER the method that is used to determine the color that a node should have
      //! @param MAX_NUMBER_LABELS the maximum number of labels that should be displayed on the dendrogram
      //! @param LABEL_MAKER the method for creating a label for nodes
      //! @param NUMBER_Y_LABELS the number of girth labels that should be displayed along the vertical axis
      //! @param SMALL_NODES_ON_BOTTOM
      //! @param MIN_MAX_GIRTH
      //! @param SCALE_WITH_NODE_SIZE
      OutputPymol
      (
        const double CYLINDER_LENGTH,
        const double CYLINDER_RADIUS,
        const double CYLINDER_SEPARATION,
        const util::ShPtr< util::FunctionInterface< Node< t_DataType, t_PrecisionType>, linal::Vector3D> > &NODE_COLORER,
        const size_t MAX_NUMBER_LABELS,
        const util::ShPtr
        <
          util::FunctionInterface
          <
            storage::Triplet //!< argument type
            <
              util::SiPtr< const Node< t_DataType, t_PrecisionType> >,
              storage::VectorND< 4, double>,
              std::string
            >,
            std::string //!< return type
          >
        > &LABEL_MAKER,
        const util::Color &COLOR,
        const size_t NUMBER_Y_LABELS,
        const bool SMALL_NODES_ON_BOTTOM,
        const storage::VectorND< 2, t_PrecisionType> &MIN_MAX_GIRTH,
        const bool SCALE_WITH_NODE_SIZE = true
      ) :
        m_CylinderLength( CYLINDER_LENGTH),
        m_CylinderRadius( CYLINDER_RADIUS),
        m_CylinderDiameter( 2.0 * CYLINDER_RADIUS),
        m_CylinderSeparation( CYLINDER_SEPARATION),
        m_NodeColorer( NODE_COLORER),
        m_MaxNodeLabels( MAX_NUMBER_LABELS),
        m_LabelMaker( LABEL_MAKER),
        m_Color( COLOR),
        m_NumberLabelsAlongYScale( NUMBER_Y_LABELS),
        m_ScaleWithNodeSize( SCALE_WITH_NODE_SIZE),
        m_SmallNodesOnBottom( SMALL_NODES_ON_BOTTOM),
        m_MinMaxGirth( MIN_MAX_GIRTH)
      {
      }

      //! @brief copy constructor
      //! @return new copy of this class
      OutputPymol< t_DataType, t_PrecisionType> *Clone() const
      {
        return new OutputPymol< t_DataType, t_PrecisionType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetCylinderLength gives "m_CylinderLength"
      //! @return returns "m_CylinderLength"
      double GetCylinderLength() const
      {
        return m_CylinderLength;
      }

      //! @brief GetCylinderRadius gives "m_CylinderRadius"
      //! @return returns "m_CylinderRadius"
      double GetCylinderRadius() const
      {
        return m_CylinderRadius;
      }

      //! @brief GetCylinderSeparation gives "m_CylinderSeparation"
      //! @return returns "m_CylinderSeparation"
      double GetCylinderSeparation() const
      {
        return m_CylinderSeparation;
      }

      //! @brief GetNodeColorer gives "m_NodeColorer"
      //! @return returns "m_NodeColorer"
      const util::ShPtr< util::FunctionInterface< Node< t_DataType, t_PrecisionType>, linal::Vector3D> > &GetNodeColorer() const
      {
        return m_NodeColorer;
      }

      //! @brief GetMaxNodeLabels gives "m_MaxNodeLabels"
      //! @return returns "m_MaxNodeLabels"
      size_t GetMaxNodeLabels() const
      {
        return m_MaxNodeLabels;
      }

      //! @brief GetLabelMaker gives "m_LabelMaker"
      //! @return returns "m_LabelMaker"
      const util::ShPtr
      <
        util::FunctionInterface
        <
          storage::Triplet //!< argument type
          <
            util::SiPtr< const Node< t_DataType, t_PrecisionType> >,
            storage::VectorND< 4, double>,
            std::string
          >,
          std::string //!< return type
        >
      > &GetLabelMaker() const
      {
        return m_LabelMaker;
      }

      //! @brief GetMaxNodeLabels gives "m_NumberLabelsAlongYScale"
      //! @return returns "m_NumberLabelsAlongYScale"
      size_t GetNumberLabelsAlongYScale() const
      {
        return m_NumberLabelsAlongYScale;
      }

      //! @brief SetCylinderLength changes "m_CylinderLength"
      //! @param NEW_CYLINDER_LENGTH the cylinder length that "m_CylinderLength" will be changed to
      void SetCylinderLength( const double NEW_CYLINDER_LENGTH)
      {
        m_CylinderLength = NEW_CYLINDER_LENGTH;
      }

      //! @brief SetCylinderRadius changes "m_CylinderRadius"
      //! @param NEW_CYLINDER_RADIUS the cylinder Radius that "m_CylinderRadius" will be changed to
      void SetCylinderRadius( const double NEW_CYLINDER_RADIUS)
      {
        m_CylinderRadius = NEW_CYLINDER_RADIUS;
      }

      //! @brief SetCylinderSeparation changes "m_CylinderSeparation"
      //! @param NEW_CYLINDER_SEPARATION the cylinder separation that "m_CylinderSeparation" will be changed to
      void SetCylinderSeparation( const double NEW_CYLINDER_SEPARATION)
      {
        m_CylinderSeparation = NEW_CYLINDER_SEPARATION;
      }

      //! @brief SetNodeColorer changes "m_NodeColorer"
      //! @param NEW_COLORER the colorer that "m_NodeColorer" will be changed to
      void SetNodeColorer
      (
        const util::ShPtr< util::FunctionInterface< Node< t_DataType, t_PrecisionType>, linal::Vector3D> > &NEW_COLORER
      )
      {
        m_NodeColorer = NEW_COLORER;
      }

      //! @brief SetMaxNodeLabels changes "m_MaxNodeLabels"
      //! @param NEW_MAX_NODE_LABELS what "m_MaxNodeLabels" will be changed to
      void SetMaxNodeLabels( const size_t NEW_MAX_NODE_LABELS)
      {
        m_MaxNodeLabels = NEW_MAX_NODE_LABELS;
      }

      //! @brief SetLabelMaker changes "m_LabelMaker"
      //! @param LABEL_MAKER what "m_LabelMaker" will be changed to
      void SetLabelMaker
      (
        const util::ShPtr
        <
          util::FunctionInterface
          <
            storage::Triplet //!< argument type
            <
              util::SiPtr< const Node< t_DataType, t_PrecisionType> >,
              storage::VectorND< 4, double>,
              std::string
            >,
            std::string //!< return type
          >
        > &LABEL_MAKER
      )
      {
        m_LabelMaker = LABEL_MAKER;
      }

      //! @brief SetNumberLabelsAlongYScale changes "m_NumberLabelsAlongYScale"
      //! @param NEW_LABELS_ALONG_Y_SCALE what "m_NumberLabelsAlongYScale" will be changed to
      void SetNumberLabelsAlongYScale( const size_t NEW_LABELS_ALONG_Y_SCALE)
      {
        m_NumberLabelsAlongYScale = NEW_LABELS_ALONG_Y_SCALE;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief WriteOutput writes this class's information about a list of nodes to a file
      //! @param FILENAME is the file to which the output will be written
      //! @param NODE_LIST the list of nodes for which the information is going to be output
      void WriteOutput
      (
        const std::string &FILENAME, const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read in members
        io::Serialize::Read( m_CylinderLength         , ISTREAM);
        io::Serialize::Read( m_CylinderRadius         , ISTREAM);
        io::Serialize::Read( m_CylinderDiameter       , ISTREAM);
        io::Serialize::Read( m_CylinderSeparation     , ISTREAM);
        io::Serialize::Read( m_NodeColorer            , ISTREAM);
        io::Serialize::Read( m_MaxNodeLabels          , ISTREAM);
        io::Serialize::Read( m_LabelMaker             , ISTREAM);
        io::Serialize::Read( m_Color                  , ISTREAM);
        io::Serialize::Read( m_NumberLabelsAlongYScale, ISTREAM);
        io::Serialize::Read( m_ScaleWithNodeSize      , ISTREAM);
        io::Serialize::Read( m_SmallNodesOnBottom     , ISTREAM);
        io::Serialize::Read( m_MinMaxGirth            , ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_CylinderLength         , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CylinderRadius         , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CylinderDiameter       , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CylinderSeparation     , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_NodeColorer            , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_MaxNodeLabels          , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_LabelMaker             , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Color                  , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_NumberLabelsAlongYScale, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ScaleWithNodeSize      , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_SmallNodesOnBottom     , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_MinMaxGirth            , OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief CalculateMinMaxGirth finds the minimum and maximum girth of a node out of a list of nodes
      //! @param "NODE_LIST" is the list of nodes which will be searched to find the maximum girth of any of the nodes
      //! @return returns a VectorND t_PrecisionType which has the minimium and maximum girths in "NODE_LIST", respectively
      storage::VectorND< 2, t_PrecisionType> CalculateMinMaxGirth
      (
        const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
      ) const;

      //! @brief DefineCylinder creates the necessary string which can be used to define a cylinder as a CGO
      //! @param START_COORDINATE start coordinates of the cylinder
      //! @param END_COORDINATE end coordinates of the cylinder
      //! @param CYLINDER_RADIUS the radius of the cylinder
      //! @param START_COLOR starting color defined by the linal::Vector3D
      //! @param END_COLOR starting color defined by the linal::Vector3D
      //! @param CGO_OBJECT_NAME the name of the dendrogram as a compiled graphics and pymol object
      //! @return returns string with Python code to create cylinder CGO object in Pymol
      static std::string DefineCylinder
      (
        const linal::Vector3D &START_COORDINATE,
        const linal::Vector3D &END_COORDINATE,
        const double CYLINDER_RADIUS,
        const linal::Vector3D &START_COLOR,
        const linal::Vector3D &END_COLOR,
        const std::string &CGO_OBJECT_NAME
      );

      //! @brief OutputNode outputs the text necessary to display a node in pymol using compiled graphics objects
      //! @param NODE the node to be displayed
      //! @param OSTREAM the ostream to which the necessary text to output "NODE" will be written
      //! @param X_TRANSLATION the position in the x direction at which "NODE" should start
      //! @param MAX_NUMBER_LABELS the number of labels that should be put onto "NODE"
      //! @param DENDROGRAM_OBJECT_NAME the name of the dendrogram object as it will be in pymol and variable in Python
      //! @param FILENAME the file name and path of the file to which the text is being written by "OSTREAM"
      //! @param BASE_GIRTH
      //! @param TOP_GIRTH
      std::ostream &OutputNode
      (
        const Node< t_DataType, t_PrecisionType> &NODE,
        std::ostream            &OSTREAM,
        double                   X_TRANSLATION,
        size_t                   MAX_NUMBER_LABELS,
        const std::string       &DENDROGRAM_OBJECT_NAME,
        const std::string       &FILENAME,
        const t_PrecisionType             BASE_GIRTH,
        const t_PrecisionType             TOP_GIRTH
      ) const;

      //! @brief DefineHeightScale creates a scale that shows the heights of the nodes in the dendrogram
      //! @param NODE_LIST the list of nodes for which the height scale will be created
      //! @param CGO_OBJECT_NAME the name of the dendrogram object as a CGO and pymol object
      //! @param OFFSET the x-position at which the scale should begin
      //! @return returns a string which has the Python code necessary to create the Y-scale as a CGO object in Pymol
      std::string DefineHeightScale
      (
        const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST, const std::string CGO_OBJECT_NAME,
        const double OFFSET
      ) const;

      //! @brief GetNodeWidth determines what the total width of a node based on the number of members it contains
      //! @param NODE_SIZE the number of members that the node contains
      //! @return a double which is the total width that a node is
      double GetNodeWidth( const size_t NODE_SIZE) const;

      //! @brief GetCylinderRadius determines the radius that a the cylinder representing a node should be
      //! @param NODE_SIZE the number of members that the node contains
      //! @return a double which is the the radius that a the cylinder representing a node should be
      double GetCylinderRadius( const size_t NODE_SIZE) const;

      //! @brief define the string for a sphere
      //! @param CENTER coordinates of the center of the sphere
      //! @param RADIUS radius
      //! @param COLOR color
      //! @param CGO_OBJECT_NAME the name of the dendrogram object as a CGO and pymol object
      std::string DefineSphere
      (
        const linal::Vector3D &CENTER,
        const double RADIUS,
        const linal::Vector3D &COLOR,
        const std::string &CGO_OBJECT_NAME
      ) const;

      size_t GetNumberBaseLeaves( const Node< t_DataType, t_PrecisionType> &NODE) const;

    private:

    }; // OutputPymol

    //! @brief WriteOutput writes this class's information about a list of nodes to a file
    //! @param FILENAME is the file to which the output will be written
    //! @param NODE_LIST the list of nodes for which the information is going to be output
    template< typename t_DataType, typename t_PrecisionType> void OutputPymol< t_DataType, t_PrecisionType>::WriteOutput
    (
      const std::string &FILENAME, const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
    ) const
    {
      // create OFStream "ofstream" - the stream to which the pymol script will be written
      io::OFStream ofstream;

      // open "ofstream" and bind it to "FILENAME"
      io::File::MustOpenOFStream( ofstream, FILENAME);

      // add the text for functions that are needed in the Python script
      ofstream
        << "from pymol.cgo import *\nfrom math import *\nfrom pymol import cmd\nfrom pymol.vfont import plain\n"
        << "cmd.bg_color(\"white\")"
        << '\n';

      // create string "dendrogram_object_name" and initialize with the name the dendrogram object will have in pymol
      // and as a variable in the Python script
      const std::string dendrogram_object_name( "dendrogram");

      // write the text to create the dendrogram as a python array
      ofstream << dendrogram_object_name << "=[]\n";

      // write the text to create the labels as a python array
      ofstream << "labels_all = []\n\n";

      // define the color for the labels
      ofstream << "# label color\n";
      ofstream << m_Color.GetName() << " = [" << m_Color->X() << ',' << m_Color->Y() << ',' << m_Color->Z() << "]\n\n";

      // create double "x_offset" and initialize with zero. This is the x-position at which the current node from
      // "NODE_LIST" will begin
      double x_offset( 0.0);

      // create size_t "number_labels_per_node" and initialize with the number of labels that should be printed for
      // each node in "NODE_LIST" so that the labels are spread evenly over all nodes in "NODE_LIST"
      const size_t number_labels_per_node( m_MaxNodeLabels / NODE_LIST.GetSize() + 1);

      // get min and max girth
      const storage::VectorND< 2, t_PrecisionType> min_max_girth( CalculateMinMaxGirth( NODE_LIST));

      // girth where all the individual leaves are
      const t_PrecisionType base_girth( m_SmallNodesOnBottom ? min_max_girth.First() : min_max_girth.Second());

      // girth where the root node is
      const t_PrecisionType top_girth( m_SmallNodesOnBottom ? min_max_girth.Second() : min_max_girth.First());

      // iterate through "NODE_LIST" in order to output each of the nodes which it contains
      for
      (
        typename util::SiPtrList< const Node< t_DataType, t_PrecisionType> >::const_iterator
          node_itr( NODE_LIST.Begin()), node_itr_end( NODE_LIST.End());
        node_itr != node_itr_end;
        ++node_itr
      )
      {
        // write the necessary text to "ofstream" to output the node currently denoted by "node_itr"
        OutputNode( **node_itr, ofstream, x_offset, number_labels_per_node, dendrogram_object_name, FILENAME, base_girth, top_girth);

        // create const size_t "node_size" and initialize with the number of members that the node currently denoted
        // by "node_itr" contains
        const size_t node_size( GetNumberBaseLeaves( **node_itr));

        // create const double "node_width" and initalize with the width of the node currently denoted by "node_itr"
        // based on "node_size", "m_CylinderDiameter", and "m_CylinderSeparation"
        const double node_width( GetNodeWidth( node_size));

        // increase "x_offset" by "node_width" so the next node will start at the correct position
        x_offset += ( node_width + m_CylinderSeparation);
      }

      // write the text to load the dendrogram object in pymol with name "dendrogram"
      ofstream << "cmd.load_cgo( " << dendrogram_object_name << ", '" << dendrogram_object_name << "')" << '\n';

      // write the text to load the labels
      ofstream << "cmd.load_cgo( labels_all, 'labels_all')\n\n";

      // write the text to show a scale in pymol giving the girths as heights in the dendrogram
      ofstream << DefineHeightScale( NODE_LIST, "y_axis", x_offset) << '\n';
    }

    //! @brief OutputNode outputs the text necessary to display a node in pymol using compiled graphics objects
    //! @param NODE the node to be displayed
    //! @param OSTREAM the ostream to which the necessary text to output "NODE" will be written
    //! @param X_TRANSLATION the position in the x direction at which "NODE" should start
    //! @param MAX_NUMBER_LABELS the number of labels that should be put onto "NODE"
    //! @param DENDROGRAM_OBJECT_NAME the name of the dendrogram object as it will be in pymol and variable in Python
    //! @param FILENAME the file name and path of the file to which the text is being written by "OSTREAM"
    //! @param BASE_GIRTH
    //! @param TOP_GIRTH
    template< typename t_DataType, typename t_PrecisionType>
    std::ostream &OutputPymol< t_DataType, t_PrecisionType>::OutputNode
    (
      const Node< t_DataType, t_PrecisionType> &NODE,
      std::ostream            &OSTREAM,
      double                   X_TRANSLATION,
      size_t                   MAX_NUMBER_LABELS,
      const std::string       &DENDROGRAM_OBJECT_NAME,
      const std::string       &FILENAME,
      const t_PrecisionType    BASE_GIRTH,
      const t_PrecisionType    TOP_GIRTH
    ) const
    {
      // create SiPtrList "expanded_node" initialize with the nodes contained in "NODE" so they are directly accessible
      const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > expanded_node( NODE.ExpandAllNodes());

      // create List of girth, node pairs to will hold the girth and the corresponding node of those in "expanded_node"
      storage::List< storage::Pair< t_PrecisionType, util::SiPtr< const Node< t_DataType, t_PrecisionType> > > > girth_node_pairs;

      // iterate through the nodes of "expanded_node" in order to fill "girth_node_map" with each girth and node
      for
      (
        typename util::SiPtrList< const Node< t_DataType, t_PrecisionType> >::const_iterator
          node_itr( expanded_node.Begin()), node_itr_end( expanded_node.End());
        node_itr != node_itr_end;
        ++node_itr
      )
      {
        // create const SiPtr "node" and initialize with the node currently denoted by "node_itr"
        const util::SiPtr< const Node< t_DataType, t_PrecisionType> > &node( *node_itr);

        // true if the girth of the node denotd by "node_itr" is defined
        if( util::IsDefined( node->GetGirth()))
        {
          // insert the girth of "node" and node into "girth_node_map"
          girth_node_pairs.PushBack
          (
            storage::Pair< t_PrecisionType, util::SiPtr< const Node< t_DataType, t_PrecisionType> > >
            (
              node->GetGirth(), node
            )
          );
        }
      }

      // create Map "node_to_x_coord_map" to hold a node and its corresponding starting x coordinate
      storage::Map< util::SiPtr< const Node< t_DataType, t_PrecisionType> >, double> node_to_x_coord_map;

      // create const size_t "total_node_size" and initialize with the number of members in "NODE"
      // This is the number of members which will be output
      const size_t total_node_size( GetNumberBaseLeaves( NODE));

      // create const double "total_node_width" and initialize with the total width that this node will be
      // based on the "total_node_size", "m_CylinderDiameter", and "m_CylinderSeparation"
      const double total_node_width( GetNodeWidth( total_node_size));

      // create const double "node_starting_x_coordinate"
      // the starting coordinate of the node is at the center of it plus the necessary "X_TRANSLATION"
      const double node_starting_x_coordinate( total_node_width / 2.0 + X_TRANSLATION);

      const linal::Vector3D rgb_color_points_base( m_NodeColorer->operator()( NODE));
      OSTREAM <<
        DefineCylinder
        (
          linal::Vector3D( node_starting_x_coordinate, NODE.GetGirth() * m_CylinderLength, 0.0), //< cylinder start
          linal::Vector3D( node_starting_x_coordinate, TOP_GIRTH * m_CylinderLength, 0.0), //< cylinder end
          GetCylinderRadius( total_node_size),
          rgb_color_points_base,
          rgb_color_points_base,
          DENDROGRAM_OBJECT_NAME
        ) << '\n';

      // insert the largest node and its position into "node_to_x_coord_map"
      node_to_x_coord_map.Insert
      (
        std::pair< util::SiPtr< const Node< t_DataType, t_PrecisionType> >, double>
        (
          girth_node_pairs.Begin()->Second(),
          node_starting_x_coordinate
        )
      );

      // create size_t "number_labels" and initialize with zero. This will keep track of the number of labels that have
      // been printed on "NODE" so that the max number of labels is not exceeded
      size_t number_labels( 0);

      // iterate through "girth_node_pairs" from the largest node to the smallest node.
      for
      (
        typename storage::List< storage::Pair< t_PrecisionType, util::SiPtr< const Node< t_DataType, t_PrecisionType> > > >::const_iterator
          node_itr( girth_node_pairs.Begin()), node_itr_end( girth_node_pairs.End());
        node_itr != node_itr_end;
        ++node_itr
      )
      {
        // create ShPtrList of Nodes "node_list" and initialize with the nodes that are in the node currently denoted
        // by "node_itr"
        const util::ShPtrList< Node< t_DataType, t_PrecisionType> > &node_list( node_itr->Second()->GetNodes());

        // create iterator "node_x_coord_itr" and initialize to the node which is denoted by "node_itr" so that the
        // x coordinate for this node can be gotten
        typename storage::Map< util::SiPtr< const Node< t_DataType, t_PrecisionType> >, double>::const_iterator node_x_coord_itr
        (
          node_to_x_coord_map.Find( node_itr->Second())
        );

        // make sure that the node could be found in "node_to_x_coord_map"
        BCL_Assert
        (
          node_x_coord_itr != node_to_x_coord_map.End(),
          "Node address not found " + util::Format()( &( *node_itr->Second())) + "\n"
          + util::Format()( node_itr->Second())
        );

        // create const double "node_x_coordinate" and initialize with the coordinate of the node denoted by
        // "node_x_coord_itr". This coordinate is used to determine the placement of the nodes which it contains.
        const double node_x_coordinate( node_x_coord_itr->second);

        // create const size_t "node_size" and initialize with the number of members in the current node
        const size_t node_size( GetNumberBaseLeaves( *node_x_coord_itr->first));

        // create const double "node_width" and initialize with the width of the current node based on "node_size"
        const double node_width( GetNodeWidth( node_size));

        // create const double "node_begin_x_coordinate" and initialize with the starting x coordinate of the current
        // node so the whole node will be centered around "node_x_coordinate"
        const double node_begin_x_coordinate( node_x_coordinate - node_width / 2.0);

        // create double "current_node_width" and initialize with zero. This will keep track of how wide the current
        // node is as its member nodes are printed.
        double current_node_width( 0.0);

        // create VectorND< 2, double> "first_last_inner_nodes_x_position" which will hold the extreme left and right
        // positions of the inner nodes. This is needed to draw the horizontal connecting line between the member nodes
        storage::VectorND< 2, double> first_last_inner_nodes_x_position
        (
          util::GetUndefined< double>(), util::GetUndefined< double>()
        );

        // create double "min_inner_node_cylinder_radius" and initialize with undefined
        // This will be used to determine the thickness of the horizontal connecting line between member nodes
        double min_inner_node_cylinder_radius( util::GetUndefined< double>());

        // iterate through the nodes (i.e. the "inner nodes") of the node denoted by "node_itr" in order to output
        // them
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::const_iterator
            inner_node_itr( node_list.Begin()), inner_node_itr_end( node_list.End());
          inner_node_itr != inner_node_itr_end;
          ++inner_node_itr, ++number_labels
        )
        {
          // create const double "inner_node_size" and initialize with the number of members in the node denoted
          // by "inner_node_itr"
          const size_t inner_node_size( GetNumberBaseLeaves( **inner_node_itr));

          // create const double "inner_node_width" and initialize with the horizontal width of the node currently
          // denoted by "inner_node_itr"
          const double inner_node_width( GetNodeWidth( inner_node_size));

          // create const double "inner_node_x_coordinate" and initialize with the x coordinate at which the current
          // node will be set i.e. where its cylinder will be placed
          const double inner_node_x_coordinate
          (
            inner_node_width / 2.0 + current_node_width + node_begin_x_coordinate
          );

          // create const double "cylinder_radius" and initialize with the radius the cylinder representing the
          // inner node will be. It is a function of the number of members contained in the inner node.
          const double cylinder_radius( GetCylinderRadius( inner_node_size));

          // true if "min_inner_node_cylinder_radius" is not defined - need to set it to "cylinder_radius"
          if( !util::IsDefined( min_inner_node_cylinder_radius))
          {
            min_inner_node_cylinder_radius = cylinder_radius;
          }

          // true if "min_inner_node_cylinder_radius" is defined and "cylinder_radius" is less
          // than "min_inner_node_cylinder_radius"
          else if( cylinder_radius < min_inner_node_cylinder_radius)
          {
            // set "min_inner_node_cylinder_radius" to "cylinder_radius"
            min_inner_node_cylinder_radius = cylinder_radius;
          }

          // create const VectorND with 3 doubles "rgb_color_points" and initialize with the color points given by
          // "m_NodeColorer"
          const linal::Vector3D rgb_color_points( m_NodeColorer->operator()( **inner_node_itr));

          // create const double "top_y_coordinate" and initialize with the top coordinate of the cylinder
          // this inner node will have: the girth of the node that contains it scaled with "m_CylinderLength"
          const double top_y_coordinate( node_itr->First());

          // create const double "bottom_y_coordinate" and initialize with the bottom coordinate of the cylinder
          // this inner node will have: its girth scaled with "m_CylinderLength"
          double bottom_y_coordinate;

          // create t_PrecisionType "node_girth" and initialize with the girth of the node denoted by "inner_node_itr"
          t_PrecisionType inner_node_girth( ( *inner_node_itr)->GetGirth());

          // true if "node_girth" is not defined - inner node is a base node
          if( inner_node_size == 1) // !util::IsDefined( inner_node_girth)
          {
            //// set "bottom_y_coordinate" to "top_y_coordinate" - this inner node won't be rendered
            //bottom_y_coordinate = top_y_coordinate;
            bottom_y_coordinate = BASE_GIRTH + math::ConvertBooleanToSign( !m_SmallNodesOnBottom) * ( math::Absolute( BASE_GIRTH - TOP_GIRTH) / ( 4.0 * m_NumberLabelsAlongYScale));
          }
          else //< "inner_node_girth" is defined
          {
            // set "bottom_y_coordinate" according to "inner_node_girth" scaled by "m_CylinderLength"
            bottom_y_coordinate = inner_node_girth;
          }

          BCL_MessageDbg
          (
            "inner_node size " + util::Format()( inner_node_size)
            + " inner width " + util::Format()( inner_node_width)
            + " inner number members " + util::Format()( ( *inner_node_itr)->GetMembers().GetSize())
            + " inner x coordinate " + util::Format()( inner_node_x_coordinate)
            + " inner number nodes " + util::Format()( ( *inner_node_itr)->GetNodes().GetSize())
            + " inner girth " + util::Format()( ( *inner_node_itr)->GetGirth())
            + " containing size " + util::Format()( node_size)
            + " containing width " + util::Format()( node_width)
            + " containing number members " + util::Format()( node_x_coord_itr->first->GetMembers().GetSize())
            + " containing x coordinate " + util::Format()( node_x_coordinate)
            + " containing number nodes " + util::Format()( node_x_coord_itr->first->GetNodes().GetSize())
            + " containing girth " + util::Format()( node_x_coord_itr->first->GetGirth())
            + " node_begin_x_coordinate " + util::Format()( node_begin_x_coordinate)
          );

          // write command to define a cylinder representing the node currently denoted by "inner_node_itr"
          OSTREAM <<
            DefineCylinder
            (
              linal::Vector3D( inner_node_x_coordinate,  bottom_y_coordinate * m_CylinderLength, 0.0), //< cylinder start
              linal::Vector3D( inner_node_x_coordinate,  top_y_coordinate * m_CylinderLength,    0.0), //< cylinder end
              cylinder_radius,
              rgb_color_points,
              rgb_color_points,
              DENDROGRAM_OBJECT_NAME
            ) << '\n';

          // true if "number_labels" is less than or equal to "MAX_NUMBER_LABELS" - need to print a label for
          // the node currently denoted by "inner_node_itr"
          if( number_labels <= MAX_NUMBER_LABELS)
          {
            // create SiPtr "inner_node" and initialize it to "inner_node_itr"
            // this it the node for which a label needs to be made
            const util::SiPtr< const Node< t_DataType, t_PrecisionType> > &inner_node( *inner_node_itr);

            // create VectorND< 4, double> "info" and inialize with data that is necessary to properly position the
            // label and make the label the correct size
            const storage::VectorND< 4, double> info
            (
              inner_node_x_coordinate, //< x coordinate of "inner_node"
              bottom_y_coordinate * m_CylinderLength, //< bottom y coordinate of "inner_node"
              top_y_coordinate * m_CylinderLength, //< top y coordinate of "inner_node"
              cylinder_radius * 2 //< diameter of "inner_node"
            );

            // make sure "m_LabelMaker" is defined
            BCL_Assert( m_LabelMaker.IsDefined(), "m_LabelMaker is not a defined pointer but should be");

            // create const std::string "label" and initialize with the text necessary to create a label
            const std::string label
            (
              m_LabelMaker->operator()
              (
                storage::Triplet
                <
                  util::SiPtr< const Node< t_DataType, t_PrecisionType> >,
                  storage::VectorND< 4, double>,
                  std::string
                >( inner_node, info, FILENAME)
              )
            );

            // write command to add a label to current inner node
            OSTREAM << label << '\n';
          }

          // insert the inner node denoted by "inner_node_itr" and its x coordinate into "node_to_x_coord_map"
          // and make sure the insert is successful
          BCL_Assert
          (
            node_to_x_coord_map.Insert
            (
              std::pair< util::SiPtr< const Node< t_DataType, t_PrecisionType> >, double>
              (
                util::ToSiPtr( **inner_node_itr), inner_node_x_coordinate
              )
            ).second,
            "Could not insert node, it must already be in map for some reason"
          );

          // increase "current_node_width" to reflect its new width based on "inner_node_width"
          current_node_width += inner_node_width + m_CylinderSeparation;

          // true if "inner_node_itr" is the first inner node to be output
          if( inner_node_itr == node_list.Begin())
          {
            // set "inner_node_x_coordinate" as the starting x coordinate of all the inner nodes of the node currently
            // denoted by "node_itr" - this is where the horizontal connecting cylinder will start
            first_last_inner_nodes_x_position.First() = inner_node_x_coordinate;
          }

          // set "inner_node_x_coordinate" as the ending x coordinate of all the inner nodes of the node currently
          // denoted by "node_itr". The last node will be the last node to set this so it will be correct at the end.
          first_last_inner_nodes_x_position.Second() = inner_node_x_coordinate;
        }

        // create VectorND of three doubles "rgb_color_points" and initialize with the color of the node denoted
        // by "node_itr". This is the color the horizontal connecting cylinder will have that connects its inner nodes
        const linal::Vector3D rgb_color_points( m_NodeColorer->operator()( *node_itr->Second()));

        // true if "first_last_inner_nodes_x_position", "first_last_inner_nodes_x_position",
        // and "min_inner_node_cylinder_radius" are all defined
        if
        (
          util::IsDefined( first_last_inner_nodes_x_position.First()) &&
          util::IsDefined( first_last_inner_nodes_x_position.Second()) &&
          util::IsDefined( min_inner_node_cylinder_radius)
        )
        {
          // true if the first and last positions to be connected by the horizontal cylinder are the same
          if( first_last_inner_nodes_x_position.First() == first_last_inner_nodes_x_position.Second())
          {
            // set the end position of the cylinder to the position of the owning node so that the inner node is
            // connected to its owning node
            first_last_inner_nodes_x_position.Second() = node_x_coordinate;
          }

          // draw the connection line between all the inner nodes of the node currently denoted by "node_itr"
          OSTREAM <<
          DefineCylinder
          (
            linal::Vector3D( first_last_inner_nodes_x_position.First(), node_itr->First() * m_CylinderLength, 0.0),
            linal::Vector3D( first_last_inner_nodes_x_position.Second(), node_itr->First() * m_CylinderLength, 0.0),
            min_inner_node_cylinder_radius,
            rgb_color_points,
            rgb_color_points,
            DENDROGRAM_OBJECT_NAME
          ) << '\n';

          // draw sphere at the end of the connecting line in order to make the output look nicer
          OSTREAM <<
          DefineSphere
          (
            linal::Vector3D( first_last_inner_nodes_x_position.First(), node_itr->First() * m_CylinderLength, 0.0),
            min_inner_node_cylinder_radius,
            rgb_color_points,
            DENDROGRAM_OBJECT_NAME
          ) << '\n';

          // draw sphere at the other end of the connecting line in order to make the output look nicer
          OSTREAM <<
          DefineSphere
          (
            linal::Vector3D( first_last_inner_nodes_x_position.Second(), node_itr->First() * m_CylinderLength, 0.0),
            min_inner_node_cylinder_radius,
            rgb_color_points,
            DENDROGRAM_OBJECT_NAME
          ) << '\n';
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief DefineCylinder creates the necessary string which can be used to define a cylinder as a CGO
    //! @param START_COORDINATE start coordinates of the cylinder
    //! @param END_COORDINATE end coordinates of the cylinder
    //! @param CYLINDER_RADIUS the radius of the cylinder
    //! @param START_COLOR starting color defined by the linal::Vector3D
    //! @param END_COLOR starting color defined by the linal::Vector3D
    //! @param CGO_OBJECT_NAME the name of the dendrogram as a compiled graphics and pymol object
    //! @return returns string with Python code to create cylinder CGO object in Pymol
    template< typename t_DataType, typename t_PrecisionType> std::string OutputPymol< t_DataType, t_PrecisionType>::DefineCylinder
    (
      const linal::Vector3D &START_COORDINATE,
      const linal::Vector3D &END_COORDINATE,
      const double CYLINDER_RADIUS,
      const linal::Vector3D &START_COLOR,
      const linal::Vector3D &END_COLOR,
      const std::string &CGO_OBJECT_NAME
    )
    {
      // initialize string stream
      std::stringstream ss_stream;

      // fill the stream up
      ss_stream
         << "cluster=[ CYLINDER, "
         << START_COORDINATE.X() << "," << START_COORDINATE.Y() << "," << START_COORDINATE.Z() << ","
         << END_COORDINATE.X() << "," << END_COORDINATE.Y() << "," << END_COORDINATE.Z() << ","
         << CYLINDER_RADIUS << ","
         << START_COLOR.X() << "," << START_COLOR.Y() << "," << START_COLOR.Z() << ","
         << END_COLOR.X() << "," << END_COLOR.Y() << "," << END_COLOR.Z()
         << "]\n" << CGO_OBJECT_NAME << ".extend( cluster)";

      // get the string from stream and return it
      return ss_stream.str();
    }

    //! @brief define the string for a sphere
    //! @param CENTER coordinates of the center of the sphere
    //! @param RADIUS radius
    //! @param COLOR color
    //! @param CGO_OBJECT_NAME the name of the dendrogram object as a CGO and pymol object
    template< typename t_DataType, typename t_PrecisionType> std::string OutputPymol< t_DataType, t_PrecisionType>::DefineSphere
    (
      const linal::Vector3D &CENTER,
      const double RADIUS,
      const linal::Vector3D &COLOR,
      const std::string &CGO_OBJECT_NAME
    ) const
    {
      // initialize string stream
      std::stringstream ss_stream;

      // fill the stream up
      ss_stream
         << "cluster=[ COLOR, " << COLOR.X() << "," << COLOR.Y() << "," << COLOR.Z()
         << ", SPHERE, " << CENTER.X() << "," << CENTER.Y() << "," << CENTER.Z() << "," << RADIUS
         << "]\n" << CGO_OBJECT_NAME << ".extend( cluster)";

      // get the string from stream and return it
      return ss_stream.str();
    }

    //! @brief DefineHeightScale creates a scale that shows the heights of the nodes in the dendrogram
    //! @param NODE_LIST the list of nodes for which the height scale will be created
    //! @param CGO_OBJECT_NAME the name of the dendrogram object as a CGO and pymol object
    //! @param OFFSET the x-position at which the scale should begin
    //! @return returns a string which has the Python code necessary to create the Y-scale as a CGO object in Pymol
    template< typename t_DataType, typename t_PrecisionType> std::string OutputPymol< t_DataType, t_PrecisionType>::DefineHeightScale
    (
      const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST, const std::string CGO_OBJECT_NAME,
      const double OFFSET
    ) const
    {
      // create t_PrecisionType "min_max_girth" and initialize with the min and max girth girth of the all the nodes
      const storage::VectorND< 2, t_PrecisionType> min_max_girth( CalculateMinMaxGirth( NODE_LIST));

      // create const t_PrecisionType "max_girth" and initialize with the max girth of all the nodes in "NODE_LIST"
      const t_PrecisionType max_girth( min_max_girth.Second());

      // create const t_PrecisionType "min_girth" and initialize with the min girth of all the nodes in "NODE_LIST"
      const t_PrecisionType min_girth( min_max_girth.First());

      // true if "min_girth" or "max_girth" is not defined - no scale will be printed
      if( !util::IsDefined( min_girth) || !util::IsDefined( max_girth) || !m_NumberLabelsAlongYScale)
      {
        // return empty string
        return std::string( "");
      }

      // create const double "max_scale_height" and initialize with the maximum height of any node in "NODE_LIST"
      const double max_scale_height( max_girth * m_CylinderLength);

      // create const double "min_scale_height" and initialize with the minimum height of any node in "NODE_LIST"
      const double min_scale_height( min_girth * m_CylinderLength);

      // create const double "text_block_size" and initialize with the size each girth label should have
      const double text_block_size( ( max_scale_height - min_scale_height) / ( m_NumberLabelsAlongYScale * 5.0));

      const double text_position_change( ( max_scale_height - min_scale_height) / m_NumberLabelsAlongYScale);

      // create const double "text_cylinder_radius" and initialize
      const double text_cylinder_radius( OutputPymolLabelString< t_PrecisionType>::GetCylinderTextRadius( text_block_size));

      // create std::string "height_scale" which will hold the text needed to create the scale as a CGO object in Pymol
      std::string height_scale;

      // add command to create the defining axes of the girth labels - determines their size and orientation
      height_scale += "axes=[[ " + util::Format()( text_block_size * 2.0) + " ,0.0,0.0],[0.0," +
        util::Format()( text_block_size) + ",0.0],[0.0,0.0,0.0]]\n";

      // add command to create Python array named "labels"
      height_scale += "labels=[]\n";

      // print out the height demarkations from the top to the bottom

      // define the color for the girth label
      height_scale += "# girth label color\n";
      height_scale += GetGirthLabelColor().GetName() + " = [";
      height_scale += util::Format()( GetGirthLabelColor()->X()) + ',';
      height_scale += util::Format()( GetGirthLabelColor()->Y()) + ',';
      height_scale += util::Format()( GetGirthLabelColor()->Z()) + "]\n\n";
      for
      (
        double label_position( max_scale_height - 0.5 * text_block_size);
        label_position >= min_scale_height - 0.5 * text_block_size;
        label_position -= text_position_change
      )
      {
        // add command to create the current girth label
        height_scale += "cyl_text(labels,plain,[" +
          util::Format()( OFFSET + 2.0 * m_CylinderDiameter + 2.0 * text_cylinder_radius) + "," +
          util::Format()( label_position) + "," + util::Format()( 0.0) + "],'" +
          util::Format().FFP( 3)( ( label_position + 0.5 * text_block_size) / m_CylinderLength) + "'," + util::Format()( text_cylinder_radius) +
          ",color=" + GetGirthLabelColor().GetName() + ",axes=axes)\n";

        // stop adding to height scale when the text position change is 0
        if( text_position_change == double( 0.0))
        {
          break;
        }
      }

      // add command to load "labels" into Pymol
      height_scale += "cmd.load_cgo(labels,'labels')\n";

      // return "height_scale" which will has the text needed to create the girth scale as a CGO object in Pymol
      return height_scale;
    }

    //! @brief CalculateMinMaxGirth finds the minimum and maximum girth of a node out of a list of nodes
    //! @param "NODE_LIST" is the list of nodes which will be searched to find the maximum girth of any of the nodes
    //! @return returns a VectorND t_PrecisionType which has the minimium and maximum girths in "NODE_LIST", respectively
    template< typename t_DataType, typename t_PrecisionType> storage::VectorND< 2, t_PrecisionType> OutputPymol< t_DataType, t_PrecisionType>::CalculateMinMaxGirth
    (
      const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
    ) const
    {
      if( util::IsDefined( m_MinMaxGirth.First()) && util::IsDefined( m_MinMaxGirth.Second()))
      {
        return m_MinMaxGirth;
      }
      else
      {
        // create t_PrecisionTypes "min_girth" and "max_girth" and initialize with undefined t_PrecisionTypes
        t_PrecisionType min_girth( util::GetUndefined< t_PrecisionType>()), max_girth( util::GetUndefined< t_PrecisionType>());

        // determine the min and max girth that is contained in "NODE_LIST"
        for
        (
          typename util::SiPtrList< const Node< t_DataType, t_PrecisionType> >::const_iterator
            node_itr( NODE_LIST.Begin()), node_itr_end( NODE_LIST.End());
          node_itr != node_itr_end;
          ++node_itr
        )
        {
          // make sure the SiPtr is defined
          BCL_Assert( node_itr->IsDefined(), "SiPtr to node is not defined but should be");

          //// create const double "current_girth" initialize with the girth of the node currently denoted by "node_itr"
          //const double current_girth( ( *node_itr)->GetGirth());

          const storage::VectorND< 2, t_PrecisionType> current_min_max( ( *node_itr)->GetSmallestLargestDefinedGirth());

          const t_PrecisionType current_min_girth( current_min_max.First());
          const t_PrecisionType current_max_girth( current_min_max.Second());

          // true if "current_girth" is defined which means we need to see if it is larger than "max_girth"
          if( util::IsDefined( current_min_girth))
          {
            // true if either ""min_girth or "max_girth" are not defined
            if( !util::IsDefined( min_girth) || !util::IsDefined( max_girth))
            {
              // set "min_girth" and "max_girth" to "current_girth"
              min_girth = current_min_girth;
              max_girth = current_min_girth;
            }

            // true if "current_girth" is smaller than "min_girth"
            else if( current_min_girth < min_girth)
            {
              // set "max_girth" to "current_girth"
              min_girth = current_min_girth;
            }
          }

          // true if "current_max_girth" is defined which means we need to see if it is larger than "max_girth"
          if( util::IsDefined( current_max_girth))
          {
            // true if either ""min_girth or "max_girth" are not defined
            if( !util::IsDefined( min_girth) || !util::IsDefined( max_girth))
            {
              // set "min_girth" and "max_girth" to "current_girth"
              min_girth = current_min_girth;
              max_girth = current_min_girth;
            }

            // true if "current_girth" is smaller than "min_girth"
            else if( current_max_girth > max_girth)
            {
              // set "max_girth" to "current_girth"
              max_girth = current_max_girth;
            }
          }
        }

        // return VectorND of t_PrecisionType initialized with the minimum and maximum girth of the all the nodes
        return storage::VectorND< 2, t_PrecisionType>( min_girth, max_girth);
      }
    }

    //! @brief GetNodeWidth determines what the total width of a node based on the number of members it contains
    //! @param NODE_SIZE the number of members that the node contains
    //! @return a double which is the total width that a node is
    template< typename t_DataType, typename t_PrecisionType> double OutputPymol< t_DataType, t_PrecisionType>::GetNodeWidth( const size_t NODE_SIZE) const
    {
      return NODE_SIZE * m_CylinderDiameter + m_CylinderSeparation * ( NODE_SIZE - 1.0);
    }

    //! @brief GetCylinderRadius determines the radius that a the cylinder representing a node should be
    //! @param NODE_SIZE the number of members that the node contains
    //! @return a double which is the the radius that a the cylinder representing a node should be
    template< typename t_DataType, typename t_PrecisionType> double OutputPymol< t_DataType, t_PrecisionType>::GetCylinderRadius( const size_t NODE_SIZE) const
    {
      if( m_ScaleWithNodeSize)
      {
        return m_CylinderRadius + 3.0 * math::Sqrt( NODE_SIZE - 1);
      }

      return m_CylinderRadius;
    }

    template< typename t_DataType, typename t_PrecisionType> size_t OutputPymol< t_DataType, t_PrecisionType>::GetNumberBaseLeaves( const Node< t_DataType, t_PrecisionType> &NODE) const
    {
      // get number of base leaves. won't be the same as number of members if the dendrogram has been cut at the bottom
      return NODE.CountNumberBaseNodes( true);
    }

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> OutputPymol< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new OutputPymol< t_DataType, t_PrecisionType>())
    );

    //! @brief operator == checks if two OutputPymols are the same
    //! @param OUTPUT_A first OutputPymols
    //! @param OUTPUT_B second OutputPymols
    //! @return returns boolean value - true if "OUTPUT_A" and "OUTPUT_B" are the same
    template< typename t_DataType, typename t_PrecisionType>
    inline bool operator ==
    (
      const OutputPymol< t_DataType, t_PrecisionType> &OUTPUT_A,
      const OutputPymol< t_DataType, t_PrecisionType> &OUTPUT_B
    )
    {
      return OUTPUT_A.GetCylinderLength()          == OUTPUT_B.GetCylinderLength()          &&
             OUTPUT_A.GetCylinderRadius()          == OUTPUT_B.GetCylinderRadius()          &&
             OUTPUT_A.GetCylinderSeparation()      == OUTPUT_B.GetCylinderSeparation()      &&
             OUTPUT_A.GetNodeColorer()             == OUTPUT_B.GetNodeColorer()             &&
             OUTPUT_A.GetMaxNodeLabels()           == OUTPUT_B.GetMaxNodeLabels()           &&
             OUTPUT_A.GetLabelMaker()              == OUTPUT_B.GetLabelMaker()              &&
             OUTPUT_A.GetNumberLabelsAlongYScale() == OUTPUT_B.GetNumberLabelsAlongYScale();
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_PYMOL_H_
