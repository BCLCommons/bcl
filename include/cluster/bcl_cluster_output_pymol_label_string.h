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

#ifndef BCL_CLUSTER_OUTPUT_PYMOL_LABEL_STRING_H_
#define BCL_CLUSTER_OUTPUT_PYMOL_LABEL_STRING_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_pymol_label_protein_model_from_string.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputPymolLabelString
    //! @brief OutputPymolLabelString is for use with the OutputPymol object and is an object that can
    //! create labels for nodes when the object being clustering is a Wrapper< string>.
    //! @details It is derived from util::FunctionInterface and takes a storage::Triplet and returns a std::string which
    //! is the text necessary to put into a Python script that can be run by Pymol in order to generate the label. The
    //! label consists of the node center, the number of members in that node, and the girth of the node.
    //!
    //! The Triplet argument to the util::FunctionIterface consists of the following :
    //! * util::SiPtr< const Node< std::string, t_PrecisionType> > is the node which will be output
    //! * storage::VectorND< 4, double> contains in this order
    //!     - the x coordinate at which the node is centered and the label will be centered
    //!     - the bottom starting y-coordinate of the node
    //!     - the top ending y-coordinate of the node
    //!     - the cylinder diameter of the node
    //! * std::string is the name of the file and path to which the the node label and all the
    //!    Python code is being written
    //!
    //! @see @link example_cluster_output_pymol_label_string.cpp @endlink
    //! @author alexanns
    //! @date September 10, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_PrecisionType>
    class OutputPymolLabelString :
      public util::FunctionInterface
      <
        storage::Triplet //!< argument type
        <
          util::SiPtr< const Node< std::string, t_PrecisionType> >,
          storage::VectorND< 4, double>,
          std::string
        >,
        std::string //!< return type
      >
    {
    protected:

    //////////
    // data //
    //////////

      //! the function used to compare the members and get the distance between members within the dendrogram
      util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
      > m_MemberDistanceFunction;

      //! the way in which the distance calculated by "m_MemberDistanceFunction" should be compared
      //! For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
      util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > m_DistanceComparator;

      //! the color for the labels
      util::Color m_Color;

      //! if true indicates that all members of the node to be labeled should be listed
      bool m_ListNodeMembers;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      OutputPymolLabelString();

      //! @brief constructor taking parameters for each member variable
      //! @param MEMBER_DISTANCE_FUNCTION the function used to compare the members and get the distance between members
      //!        within the dendrogram
      //! @param DISTANCE_COMPARATOR the way in which the distance calculated by "m_MemberDistanceFunction" should
      //!        be compared. For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
      //! @param COLOR the color to use for the labels
      //! @param LIST_MEMBERS if true indicates that all members of the node to be labeled should be listed
      OutputPymolLabelString
      (
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
        > &MEMBER_DISTANCE_FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &DISTANCE_COMPARATOR,
        const util::Color &COLOR,
        const bool LIST_MEMBERS = false
      );

      //! @brief virtual copy constructor
      //! @return new copy of this class
      OutputPymolLabelString *Clone() const;

      //! @brief virtual destructor
      virtual ~OutputPymolLabelString();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking a Triplet and returns a std::string containing the label text
      //! @param TRIPLET the Triplet to the util::FunctionIterface consists of the following :
      //! * util::SiPtr< const Node< std::string, t_PrecisionType> > is the node which will be output
      //! * storage::VectorND< 4, double> contains in this order
      //!     - the x coordinate at which the node is centered and the label will be centered
      //!     - the bottom starting y-coordinate of the node
      //!     - the top ending y-coordinate of the node
      //!     - the cylinder diameter of the node
      //! * std::string is the name of the file and path to which the the node label and all the
      //!    Python code is being written
      //! @return returns a string which is the Python code necessary to create the label in Pymol
      std::string operator()
      (
        const storage::Triplet
        <
          util::SiPtr< const Node< std::string, t_PrecisionType> >,
          storage::VectorND< 4, double>,
          std::string
        > &TRIPLET
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief GetCylinderText creates text to create a text compiled graphics object cylinder in Pymol using Python
      //! @param TEXT_OBJECT_NAME the name of the text object in the Python code and in Pymol
      //! @param X_COORDINATE the x coordinate where the text will be placed
      //! @param Y_COORDINATE the y coordinate where the text will be placed
      //! @param Z_COORDINATE the z coordinate where the text will be placed
      //! @param LABEL_TEXT the actual text that will show up in pymol to label the node
      //! @param TEXT_CYLINDER_RADIUS the radius of the cylinders that will make up the text
      //! @param TEXT_BLOCK_SIZE the total block size of the text
      //! @param COLOR the red,green,blue color point the text will have
      //! @return returns a string which can create a cylinder text object in Pymol using Python
      static std::string GetCylinderText
      (
        const std::string &TEXT_OBJECT_NAME,
        const double X_COORDINATE,
        const double Y_COORDINATE,
        const double Z_COORDINATE,
        const std::string &LABEL_TEXT,
        const double TEXT_CYLINDER_RADIUS,
        const double TEXT_BLOCK_SIZE,
        const util::Color &COLOR
      );

      //! @brief GetLabelCommand creates a the Python code to make a basic CGO text label for a Node in Pymol
      //!        The label consists of the size (i.e. number of members) and girth of the node
      //!        The labels are placed starting at the point 1/4 of the way between the provided y-coordinates below
      //!        one another
      //! @param NODE the node for which the label will be made
      //! @param LABEL_OBJECT_NAME the name the CGO object will have in Pymol and as a variable in Python
      //! @param X_COORDINATE the x coordinate where the labels will be placed
      //! @param Y_COORDINATE_TOP the top y-coordinate of the node
      //! @param Y_COORDINATE_BOTTOM the bottom y-coordinate of the node
      //! @param NODE_CYLINDER_DIAMETER the diameter of the cylinder representing "NODE"
      //! @param COLOR the red,green,blue color point the text will have
      //! @return returns the text which is the Python code to make a basic CGO text label for a Node in Pymol
      template< typename t_DataType>
      static std::string GetLabelCommand
      (
        const Node< t_DataType, t_PrecisionType> &NODE, const std::string &LABEL_OBJECT_NAME, const double X_COORDINATE,
        const double Y_COORDINATE_TOP, const double Y_COORDINATE_BOTTOM, const double NODE_CYLINDER_DIAMETER,
        const util::Color &COLOR
      );

      //! @brief GetTextZCoordinate gives the z coordinate where the text for a label should start
      //!        Since the dendrogram is centered on z=0, only the diameter of the node and the size of the text
      //!        is needed to determine where the text label should start
      //! @param NODE_CYLINDER_DIAMETER the diameter of the node for which the label will be applied
      //! @param TEXT_CYLINDER_RADIUS the radial thickness of the cylinders that will make up the text
      //! @return returns a double which is the z coordinate where the text should start
      static double GetTextZCoordinate( const double NODE_CYLINDER_DIAMETER, const double TEXT_CYLINDER_RADIUS);

      //! @brief GetTextBlockSize gives the size that a piece of text should have
      //! @param Y_COORDINATE_TOP the top y-coordinate of the node
      //! @param Y_COORDINATE_BOTTOM the bottom y-coordinate of the node
      //! @return returns a double which is the size the text block should be in order to fit on its node's cylinder
      static double GetTextBlockSize( const double Y_COORDINATE_TOP, const double Y_COORDINATE_BOTTOM);

      //! @brief GetCylinderTextRadius gives the radius of the cylinders that create a CGO text
      //! @param TEXT_BLOCK_SIZE the total size of the text block
      //! @return returns a double which is the radius of the cylinders that create a CGO text
      static double GetCylinderTextRadius( const double TEXT_BLOCK_SIZE);

    private:

      //! @brief GetUniqueObjectName creates a name for a label object for a node that will be unique to other nodes
      //!        It will not be unique for multiple labels of the same node; a prefix or postfix must be added to the
      //!        returned value to accomplish this
      //! @param X_COORDINATE the x coordinate of the node for which the label will be made
      //! @param GIRTH the girth of the node for which the label will be made
      //! @return returns a std::string which is a combination of "X_COORDINATE" and "GIRTH"
      //!         Since no two nodes should have the same x coordinate and girth, this string should be unique
      static std::string GetUniqueObjectName( const double X_COORDINATE, const double GIRTH);

    }; // OutputPymolLabelString

    //! @brief GetLabelCommand creates a the Python code to make a basic CGO text label for a Node in Pymol
    //!        The label consists of the size (i.e. number of members) and girth of the node
    //!        The labels are placed starting at the point 1/4 of the way between the provided y-coordinates below
    //!        one another
    //! @param NODE the node for which the label will be made
    //! @param LABEL_OBJECT_NAME the name the CGO object will have in Pymol and as a variable in Python
    //! @param X_COORDINATE the x coordinate where the labels will be placed
    //! @param Y_COORDINATE_TOP the top y-coordinate of the node
    //! @param Y_COORDINATE_BOTTOM the bottom y-coordinate of the node
    //! @param NODE_CYLINDER_DIAMETER the diameter of the cylinder representing "NODE"
    //! @param COLOR the red,green,blue color point the text will have
    //! @return returns the text which is the Python code to make a basic CGO text label for a Node in Pymol
    template< typename t_PrecisionType>
    template< typename t_DataType> std::string OutputPymolLabelString< t_PrecisionType>::GetLabelCommand
    (
      const Node< t_DataType, t_PrecisionType> &NODE, const std::string &LABEL_OBJECT_NAME, const double X_COORDINATE,
      const double Y_COORDINATE_TOP, const double Y_COORDINATE_BOTTOM, const double NODE_CYLINDER_DIAMETER, const util::Color &COLOR
    )
    {
      // create std::string "label" which will hold the necessary text to create the commands
      std::string label( "");

      // create const double "text_size" which is the total height of the text block and initialize with the text
      // block being one twentieth of the total length of the node
      const double text_size( GetTextBlockSize( Y_COORDINATE_TOP, Y_COORDINATE_BOTTOM));

      // create const double "text_cylinder_radius" which holds the radius of the cylinders that create the text
      // and initialize with one eighth of "text_size"
      const double text_cylinder_radius( GetCylinderTextRadius( text_size));

      // add command for creating the Python array to hold the label to "label"
      label += LABEL_OBJECT_NAME + "=[]\n";

      // create const double "midpoint" and initialize with the halfway point of "node" in the y-direction
      const double midpoint( ( Y_COORDINATE_BOTTOM + Y_COORDINATE_TOP) / 2.0);

      // create const double "quarter_point" and initialize with the point halfway between the bottom of "node"
      // and the "midpoint" of "node". This is where the label text will be positioned.
      const double quarter_point( ( midpoint + Y_COORDINATE_BOTTOM) / 2.0);

      // create const double "z_coordinate" and initialize with the starting z coordinate of the text
      const double z_coordinate( GetTextZCoordinate( NODE_CYLINDER_DIAMETER, text_cylinder_radius));

      // add command to print out the size of the cluster to "label"
      label += GetCylinderText
        (
          LABEL_OBJECT_NAME, X_COORDINATE, quarter_point, z_coordinate,
          util::Format()( NODE.GetMembers().GetSize()), text_cylinder_radius, text_size, COLOR
        );

      // create double "girth" and initialize with the girth of "NODE"
      std::string girth_str( util::Format()( NODE.GetGirth()));

      // true if the girth is not defined
      if( !util::IsDefined( NODE.GetGirth()))
      {
        // don't want to work with undefined numbers so set the girth to zero
        girth_str = "undef";
      }

      // add command to print out the girth of the cluster to "label"
      label += GetCylinderText
      (
        LABEL_OBJECT_NAME, X_COORDINATE, quarter_point - ( text_size + 2 * text_cylinder_radius),
        z_coordinate, girth_str, text_cylinder_radius, text_size, COLOR
      );

      // add command to print out the node identifier
      label += GetCylinderText
      (
        LABEL_OBJECT_NAME, X_COORDINATE, quarter_point + ( text_size + 2 * text_cylinder_radius),
        z_coordinate, util::Format()( NODE.GetIdentifier()), text_cylinder_radius, text_size, COLOR
      );

      // add command to load the label object in Pymol
      label += "labels_all.extend( " + LABEL_OBJECT_NAME + ")\n";

      // return "label" which has all the text needed to create the label for "NODE"
      return label;
    }

    //! single instance of that class
    template< typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> OutputPymolLabelString< t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new OutputPymolLabelString< t_PrecisionType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_PrecisionType>
    OutputPymolLabelString< t_PrecisionType>::OutputPymolLabelString() :
      m_MemberDistanceFunction(),
      m_DistanceComparator(),
      m_Color( util::GetColors().e_Magenta),
      m_ListNodeMembers( false)
    {
    }

    //! @brief constructor taking parameters for each member variable
    //! @param MEMBER_DISTANCE_FUNCTION the function used to compare the members and get the distance between members
    //!        within the dendrogram
    //! @param DISTANCE_COMPARATOR the way in which the distance calculated by "m_MemberDistanceFunction" should
    //!        be compared. For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
    //! @param COLOR the color to use for the labels
    //! @param LIST_MEMBERS if true indicates that all members of the node to be labeled should be listed
    template< typename t_PrecisionType>
    OutputPymolLabelString< t_PrecisionType>::OutputPymolLabelString
    (
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
      > &MEMBER_DISTANCE_FUNCTION,
      const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &DISTANCE_COMPARATOR,
      const util::Color &COLOR,
      const bool LIST_MEMBERS
    ) :
      m_MemberDistanceFunction( MEMBER_DISTANCE_FUNCTION),
      m_DistanceComparator( DISTANCE_COMPARATOR),
      m_Color( COLOR),
      m_ListNodeMembers( LIST_MEMBERS)
    {
    }

    //! @brief virtual copy constructor
    //! @return new copy of this class
    template< typename t_PrecisionType>
    OutputPymolLabelString< t_PrecisionType> *OutputPymolLabelString< t_PrecisionType>::Clone() const
    {
      return new OutputPymolLabelString< t_PrecisionType>( *this);
    }

    //! @brief virtual destructor
    template< typename t_PrecisionType>
    OutputPymolLabelString< t_PrecisionType>::~OutputPymolLabelString()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_PrecisionType>
    const std::string &OutputPymolLabelString< t_PrecisionType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking a Triplet and returns a std::string containing the label text
    //! @param TRIPLET the Triplet to the util::FunctionIterface consists of the following :
    //! * util::SiPtr< const Node< std::string, t_PrecisionType> > is the node which will be output
    //! * storage::VectorND< 4, double> contains in this order
    //!     - the x coordinate at which the node is centered and the label will be centered
    //!     - the bottom starting y-coordinate of the node
    //!     - the top ending y-coordinate of the node
    //!     - the cylinder diameter of the node
    //! * std::string is the name of the file and path to which the the node label and all the
    //!    Python code is being written
    //! @return returns a string which is the Python code necessary to create the label in Pymol
    template< typename t_PrecisionType>
    std::string OutputPymolLabelString< t_PrecisionType>::operator()
    (
      const storage::Triplet
      <
        util::SiPtr< const Node< std::string, t_PrecisionType> >,
        storage::VectorND< 4, double>,
        std::string
      > &TRIPLET
    ) const
    {
      // create util::SiPtr< const Node< std::string, t_PrecisionType> > "node" and initialize from "TRIPLET"
      const util::SiPtr< const Node< std::string, t_PrecisionType> > node( TRIPLET.First());

      // create double "x_coordinate" and initialize from "TRIPLET"
      const double x_coordinate( TRIPLET.Second()( 0));

      // create double "y_coordinate_bottom" and initialize from "TRIPLET"
      const double y_coordinate_bottom( TRIPLET.Second()( 1));

      // create double "y_coordinate_top" and initialize from "TRIPLET"
      const double y_coordinate_top( TRIPLET.Second()( 2));

      // create double "y_coordinate_top" and initialize from "TRIPLET"
      const double node_cylinder_diameter( TRIPLET.Second()( 3));

      // create util::SiPtr< const std::string> "script_filename" and initialize from "TRIPLET"
      const std::string script_filename( TRIPLET.Third());

      // create SiPtr to std::string"node_center" and initialize with the node member that
      // is closest to all the other members in "node"
      const util::SiPtr< const std::string> node_center
      (
        node->GetCenter( m_MemberDistanceFunction, m_DistanceComparator)
      );

      // create const double "girth" and initialize with the girth of "node"
      double girth( node->GetGirth());

      // create string "label" which will hold the text necessary to create the entire label for "node"
      // initialize as empty string
      std::string label( "");

      // create string "object_name" and initialize with the x coordinate and girth
      // this is the name this label will have when loaded into Pymol
      const std::string object_name( GetUniqueObjectName( x_coordinate, girth));

      // add command to print the girth and size of "node"
      label += GetLabelCommand< std::string>
      (
        *node, "size_girth" + object_name, x_coordinate,
        y_coordinate_top, y_coordinate_bottom, node_cylinder_diameter, m_Color
      );

      // create const std::string "file_object_name" and initialize with the name the node center member string
      // object will have when loaded into Pymol
      const std::string file_object_name( "center_" + object_name);

      // add command to create Python object named "file_object_name"
      label += file_object_name + "=[]\n";

      // create const double "text_block_size" and initialize with the results of GetTextBlockSize
      const double text_block_size( GetTextBlockSize( y_coordinate_top, y_coordinate_bottom));

      // create const double "text_cylinder_radius" and initialize with the radius of the cylinders that will make
      // up the text that will be the node center member object string
      const double text_cylinder_radius( GetCylinderTextRadius( text_block_size));

      // add command to get a Pymol object which is the text of the string that is the center member object of "node"
      label += GetCylinderText
      (
        file_object_name,
        x_coordinate,
        ( ( ( y_coordinate_top + y_coordinate_bottom) / 2) + y_coordinate_top) / 2,
        GetTextZCoordinate( node_cylinder_diameter, text_cylinder_radius), *node_center,
        text_cylinder_radius,
        text_block_size, m_Color
      );

      // add command to load the label object in Pymol
      label += "labels_all.extend( " + file_object_name + ")\n";

      if( m_ListNodeMembers)
      {
        OutputPymolLabelProteinModelFromString< t_PrecisionType>::ListNodeMembers
        (
          *node,
          io::File::RemoveFullExtension( script_filename) + "_node_" + util::Format()( node->GetIdentifier()) + ".ls"
        );
      }

      // return "label"
      return label;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_PrecisionType>
    std::istream &OutputPymolLabelString< t_PrecisionType>::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MemberDistanceFunction, ISTREAM);
      io::Serialize::Read( m_DistanceComparator    , ISTREAM);
      io::Serialize::Read( m_Color                 , ISTREAM);
      io::Serialize::Read( m_ListNodeMembers       , ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    template< typename t_PrecisionType>
    std::ostream &OutputPymolLabelString< t_PrecisionType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MemberDistanceFunction, OSTREAM, INDENT);
      io::Serialize::Write( m_DistanceComparator    , OSTREAM, INDENT);
      io::Serialize::Write( m_Color                 , OSTREAM, INDENT);
      io::Serialize::Write( m_ListNodeMembers       , OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetCylinderText creates text necessary to create a text compiled graphics object in Pymol using Python
    //! @param TEXT_OBJECT_NAME the name of the text object in the Python code and in Pymol
    //! @param X_COORDINATE the x coordinate where the text will be placed
    //! @param Y_COORDINATE the y coordinate where the text will be placed
    //! @param Z_COORDINATE the z coordinate where the text will be placed
    //! @param LABEL_TEXT the actual text that will show up in pymol to label the node
    //! @param TEXT_CYLINDER_RADIUS the radius of the cylinders that will make up the text
    //! @param TEXT_BLOCK_SIZE the total block size of the text
    //! @param COLOR the red,green,blue color point the text will have
    //! @return returns a string which can create a cylinder text object in Pymol using Python
    template< typename t_PrecisionType>
    std::string OutputPymolLabelString< t_PrecisionType>::GetCylinderText
    (
      const std::string &TEXT_OBJECT_NAME,
      const double X_COORDINATE,
      const double Y_COORDINATE,
      const double Z_COORDINATE,
      const std::string &LABEL_TEXT,
      const double TEXT_CYLINDER_RADIUS,
      const double TEXT_BLOCK_SIZE,
      const util::Color &COLOR
    )
    {
      // create string "text" which will hold the text necessary to write out "LABEL_TEXT" according to the
      // parameters passed. Initialize as an empty string.
      std::string text( "");

      // add to "text" the command determining the orientation of the text block
      text +=
        "axes=[[  0.0,0.0," + util::Format()( TEXT_BLOCK_SIZE) + "],[0.0,"
        + util::Format()( TEXT_BLOCK_SIZE) + ",0.0],[0.0,0.0,0.0]]\n";

      // add to "text" the command for printing the text
      text +=
        "cyl_text(" + TEXT_OBJECT_NAME + ",plain,["
        + util::Format()( X_COORDINATE) //< x coordinate of label
        + "," + util::Format()( Y_COORDINATE) //< y coordinate of label
        + "," + util::Format()( Z_COORDINATE) + "],'" //< z coordinate of the label
        + LABEL_TEXT + "'," //< print the actual text
        + util::Format()( TEXT_CYLINDER_RADIUS) //< set the radius of the cylinders making up the text
        + ",color=" + COLOR.GetName()
        + ",axes=axes)\n";

      // return "text"
      return text;
    }

    //! @brief GetTextZCoordinate gives the z coordinate where the text for a label should start
    //!        Since the dendrogram is centered on z=0, only the diameter of the node and the size of the text
    //!        is needed to determine where the text label should start
    //! @param NODE_CYLINDER_DIAMETER the diameter of the node for which the label will be applied
    //! @param TEXT_CYLINDER_RADIUS the radial thickness of the cylinders that will make up the text
    //! @return returns a double which is the z coordinate where the text should start
    template< typename t_PrecisionType>
    double OutputPymolLabelString< t_PrecisionType>::GetTextZCoordinate
    (
      const double NODE_CYLINDER_DIAMETER, const double TEXT_CYLINDER_RADIUS
    )
    {
      // return the z coordinate where the text should start so it is outside the dendrogram
      return NODE_CYLINDER_DIAMETER / 1.75 + 2.0 * TEXT_CYLINDER_RADIUS;
    }

    //! @brief GetTextBlockSize gives the size that a piece of text should have
    //! @param Y_COORDINATE_TOP the top y-coordinate of the node
    //! @param Y_COORDINATE_BOTTOM the bottom y-coordinate of the node
    //! @return returns a double which is the size the text block should be in order to fit on its node's cylinder
    template< typename t_PrecisionType>
    double OutputPymolLabelString< t_PrecisionType>::GetTextBlockSize( const double Y_COORDINATE_TOP, const double Y_COORDINATE_BOTTOM)
    {
      // return the size the text should have so it fits with its node's cylinder
      return std::abs( Y_COORDINATE_TOP - Y_COORDINATE_BOTTOM) / double( 20.0);
    }

    //! @brief GetCylinderTextRadius gives the radius of the cylinders that create a CGO text
    //! @param TEXT_BLOCK_SIZE the total size of the text block
    //! @return returns a double which is the radius of the cylinders that create a CGO text
    template< typename t_PrecisionType>
    double OutputPymolLabelString< t_PrecisionType>::GetCylinderTextRadius( const double TEXT_BLOCK_SIZE)
    {
      return TEXT_BLOCK_SIZE / double( 8.0);
    }

    //! @brief GetUniqueObjectName creates a name for a label object for a node that will be unique to other nodes
    //!        It will not be unique for multiple labels of the same node; a prefix or postfix must be added to the
    //!        returned value to accomplish this
    //! @param X_COORDINATE the x coordinate of the node for which the label will be made
    //! @param GIRTH the girth of the node for which the label will be made
    //! @return returns a std::string which is a combination of "X_COORDINATE" and "GIRTH"
    //!         Since no two nodes should have the same x coordinate and girth, this string should be unique
    template< typename t_PrecisionType>
    std::string OutputPymolLabelString< t_PrecisionType>::GetUniqueObjectName( const double X_COORDINATE, const double GIRTH)
    {
      // create double "girth" and initialize with "GIRTH"
      double girth( GIRTH);

      // true if "GIRTH" is undefined - need to set "girth" to zero so that strings for undefined numbers
      // aren't used in the file name
      if( !util::IsDefined( GIRTH))
      {
        // set girth to zero
        girth = 0.0;
      }

      // return string which can be the object name for a label loaded in Pymol or a Python variable
      return "_" + util::ReplaceString( util::Format()( X_COORDINATE), ".", "_") + "_" +
        util::ReplaceString( util::Format()( girth), ".", "_");
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_PYMOL_LABEL_STRING_H_ 
