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

#ifndef BCL_CLUSTER_OUTPUT_PYMOL_LABEL_PROTEIN_MODEL_FROM_STRING_H_
#define BCL_CLUSTER_OUTPUT_PYMOL_LABEL_PROTEIN_MODEL_FROM_STRING_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_pymol_label_small_molecule.h"
#include "bcl_cluster_output_pymol_label_string.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_line.h"
#include "pdb/bcl_pdb_residue.h"
#include "score/bcl_score_radius_of_gyration.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputPymolLabelProteinModelFromString
    //! @brief OutputPymolLabelProteinModelFromString is for use with the OutputPymol object and is an object that can
    //! create labels for nodes when the object being clustering is a string that denotes a PDB file for a  protein
    //! model
    //! @details It is derived from util::FunctionInterface and takes a storage::Triplet and returns a std::string which
    //! is the text necessary to put into a Python script that can be run by Pymol in order to generate the label
    //! The label consists of the
    //! actual protein which is the node center, the number of members in that node, and the girth of the node
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
    //! @see @link example_cluster_output_pymol_label_protein_model_from_string.cpp @endlink
    //! @author alexanns
    //! @date September 3, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_PrecisionType>
    class OutputPymolLabelProteinModelFromString :
      public util::FunctionInterface
      <
        storage::Triplet // argument type
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

      //! the prefix which must be added to the clustered Wrapper< string> member in order to create the true
      //! PDB filename
      std::string m_Prefix;

      //! the prefix which must be added to the clustered Wrapper< string> member in order to create the true
      //! PDB filename
      std::string m_Postfix;

      //! lable color
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
      OutputPymolLabelProteinModelFromString();

      //! @brief constructor taking parameters for each member variable
      //! @param MEMBER_DISTANCE_FUNCTION the function used to compare the members and get the distance between members
      //!        within the dendrogram
      //! @param DISTANCE_COMPARATOR the way in which the distance calculated by "m_MemberDistanceFunction" should
      //!        be compared. For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
      //! @param PREFIX the prefix which must be added to the clustered Wrapper< string> member in order to create the
      //!        true PDB filename
      //! @param POSTFIX the postfix which must be added to the clustered Wrapper< string> member in order to create
      //!        the true PDB filename
      //! @param COLOR label color
      //! @param LIST_MEMBERS if true indicates that all members of the node to be labeled should be listed
      OutputPymolLabelProteinModelFromString
      (
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
        > &MEMBER_DISTANCE_FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &DISTANCE_COMPARATOR,
        const std::string &PREFIX,
        const std::string &POSTFIX,
        const std::string &COLOR,
        const bool LIST_MEMBERS = false
      );

      //! @brief virtual copy constructor
      //! @return new copy of this class
      OutputPymolLabelProteinModelFromString *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief GetProteinModelFromString creates a protein model from a filename
      //! @param PDB_FILENAME the name of the PDB file from which the protein model will be created
      //! @return returns a protein model that has been created from "PDB_FILENAME"
      static assemble::ProteinModel GetProteinModelFromString
      (
        const std::string &PDB_FILENAME
      );

      //! @brief writes all the members of the given node to individual pdb files
      //! @param NODE the node whose members will be written into pdb files
      //! @param LIST_FILENAME which will contain the list of pdbs that are members of the given node
      //! @param PREFIX the prefix that will be added to the member name
      //! @param POSTFIX the postfix that will be added to the member name
      static void ListNodeMembers
      (
        const Node< std::string, t_PrecisionType> &NODE,
        const std::string &LIST_FILENAME,
        const std::string &PREFIX = "",
        const std::string &POSTFIX = ""
      );

      //! @brief removes characters from string which could cause problems within python code
      //! @param STRING the string which will be cleaned
      static void CleanString( std::string &STRING);

    }; // OutputPymolLabelProteinModelFromString

    //! single instance of that class
    template< typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> OutputPymolLabelProteinModelFromString< t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new OutputPymolLabelProteinModelFromString< t_PrecisionType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_PrecisionType>
    OutputPymolLabelProteinModelFromString< t_PrecisionType>::OutputPymolLabelProteinModelFromString() :
      m_MemberDistanceFunction(),
      m_DistanceComparator(),
      m_Prefix(),
      m_Postfix(),
      m_Color( util::GetColors().e_Magenta),
      m_ListNodeMembers( false)
    {
    }

    //! @brief constructor taking parameters for each member variable
    //! @param MEMBER_DISTANCE_FUNCTION the function used to compare the members and get the distance between members
    //!        within the dendrogram
    //! @param DISTANCE_COMPARATOR the way in which the distance calculated by "m_MemberDistanceFunction" should
    //!        be compared. For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
    //! @param PREFIX the prefix which must be added to the clustered Wrapper< string> member in order to create the
    //!        true PDB filename
    //! @param POSTFIX the postfix which must be added to the clustered Wrapper< string> member in order to create
    //!        the true PDB filename
    //! @param COLOR label color
    //! @param LIST_MEMBERS if true indicates that all members of the node to be labeled should be listed
    template< typename t_PrecisionType>
    OutputPymolLabelProteinModelFromString< t_PrecisionType>::OutputPymolLabelProteinModelFromString
    (
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
      > &MEMBER_DISTANCE_FUNCTION,
      const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &DISTANCE_COMPARATOR,
      const std::string &PREFIX,
      const std::string &POSTFIX,
      const std::string &COLOR,
      const bool LIST_MEMBERS
    ) :
      m_MemberDistanceFunction( MEMBER_DISTANCE_FUNCTION),
      m_DistanceComparator( DISTANCE_COMPARATOR),
      m_Prefix( PREFIX),
      m_Postfix( POSTFIX),
      m_Color( COLOR),
      m_ListNodeMembers( LIST_MEMBERS)
    {
    }

    //! @brief virtual copy constructor
    //! @return new copy of this class
    template< typename t_PrecisionType>
    OutputPymolLabelProteinModelFromString< t_PrecisionType> *OutputPymolLabelProteinModelFromString< t_PrecisionType>::Clone() const
    {
      return new OutputPymolLabelProteinModelFromString< t_PrecisionType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_PrecisionType>
    const std::string &OutputPymolLabelProteinModelFromString< t_PrecisionType>::GetClassIdentifier() const
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
    std::string OutputPymolLabelProteinModelFromString< t_PrecisionType>::operator()
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
      const std::string &script_filename( TRIPLET.Third());

      // create SiPtr to Wrapper string "node_center" and initialize with the protein model that is closest to the
      // all the other protein models in the node
      const util::SiPtr< const std::string> node_center
      (
        node->GetCenter( m_MemberDistanceFunction, m_DistanceComparator)
      );

      // create string "pdb_filename" and initialize with the name of the pdb file that the moved protein model will
      // be written to (moved to correct position in dendrogram)
      // "pdb_filename" also contains a path for the PDB file which is the same as the path of "script_filename"
      const std::string pdb_filename( m_Prefix + *node_center + m_Postfix);

      // create string "label" which will hold the text necessary to create the entire label for "node"
      std::string label;

      // create object name for protein in pymol
      std::string pymol_obj_name( io::File::RemoveFullExtension( *node_center) + "_" + util::Format()( node->GetIdentifier()));

      CleanString( pymol_obj_name);

      // load the pdb file into pymol
      label += OutputPymolLabelSmallMolecule< t_PrecisionType>::LoadPymolObject( pdb_filename, pymol_obj_name);

      // create const double "sqr_radius_gyration" and initialize with the square radius of gyration of "model_copy"
      const double sqr_radius_gyration( math::Sqrt( score::RadiusOfGyration().SquareRadiusOfGyration( GetProteinModelFromString( pdb_filename))));

      // move the molecule to the position of the node
      const double translate_x_coordinate( x_coordinate - 2.0 * sqr_radius_gyration - node_cylinder_diameter / 2.0);
      const double translate_y_coordinate( ( y_coordinate_top + y_coordinate_bottom) / 2.0 - 2.0 * sqr_radius_gyration);
      const double translate_z_coordinate( node_cylinder_diameter / 2.0 + 2.0 * sqr_radius_gyration);
      label += OutputPymolLabelSmallMolecule< t_PrecisionType>::TranslatePymolObject
        ( translate_x_coordinate, translate_y_coordinate, translate_z_coordinate, pymol_obj_name);

      // add command to show the protein as cartoon "label"
      label += "cmd.show(\"cartoon\", \"" + pymol_obj_name + "\" )\n";

      // add command to hide the protein sticks "label"
      label += "cmd.hide(\"lines\", \"" + pymol_obj_name + "\" )\n";

      // add command to show the protein with each chain colored rainbow to "label"
      label += "util.chainbow(\"" + pymol_obj_name + "\")\n";

      // add command to show the size and girth of "node"
      label += OutputPymolLabelString< t_PrecisionType>::template GetLabelCommand< std::string>
      (
        *node, "size_girth_current_label",
        x_coordinate, y_coordinate_top, y_coordinate_bottom, node_cylinder_diameter, m_Color
      );

      // label the node with the pdb name of the pdb which is the center of the node
      label += OutputPymolLabelSmallMolecule< t_PrecisionType>::GetObjectNameLabel
      (
        x_coordinate, y_coordinate_top, y_coordinate_bottom, node_center, node_cylinder_diameter, m_Color
      );

      if( m_ListNodeMembers)
      {
        ListNodeMembers
        (
          *node,
          io::File::RemoveFullExtension( script_filename) + "_node_" + util::Format()( node->GetIdentifier()) + ".ls",
          m_Prefix,
          m_Postfix
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
    std::istream &OutputPymolLabelProteinModelFromString< t_PrecisionType>::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MemberDistanceFunction, ISTREAM);
      io::Serialize::Read( m_DistanceComparator    , ISTREAM);
      io::Serialize::Read( m_Prefix                , ISTREAM);
      io::Serialize::Read( m_Postfix               , ISTREAM);
      io::Serialize::Read( m_Color                 , ISTREAM);
      io::Serialize::Read( m_ListNodeMembers       , ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    template< typename t_PrecisionType>
    std::ostream &OutputPymolLabelProteinModelFromString< t_PrecisionType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MemberDistanceFunction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceComparator    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Prefix                , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Postfix               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Color                 , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ListNodeMembers       , OSTREAM, INDENT) << '\n';

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetProteinModelFromString creates a protein model from a filename
    //! @param PDB_FILENAME the name of the PDB file from which the protein model will be created
    //! @return returns a protein model that has been created from "PDB_FILENAME"
    template< typename t_PrecisionType>
    assemble::ProteinModel OutputPymolLabelProteinModelFromString< t_PrecisionType>::GetProteinModelFromString
    (
      const std::string &PDB_FILENAME
    )
    {
      // create "factory" to create protein model with amino acids of type AABackBone
      const pdb::Factory factory;

      // initialize write and read stream objects
      io::IFStream read;

      // open pdb file
      io::File::MustOpenIFStream( read, PDB_FILENAME);

      // true is used to advise handler to ignore clashes in the structure and insert residues as they are suggested
      // by the numbering without regard for the bond information
      pdb::Handler pdb( read, true);

      // close "read"
      io::File::CloseClearFStream( read);

      // return a protein model which is created by "factory" based on "pdb"
      return factory.ProteinModelFromPDB( pdb);
    }

    //! @brief writes all the members of the given node to individual mdl files
    //! @param NODE the node whose members will be written into mdl files
    //! @param LIST_FILENAME which will contain the list of pdbs that are members of the given node
    //! @param PREFIX the prefix that will be added to the member name
    //! @param POSTFIX the postfix that will be added to the member name
    template< typename t_PrecisionType>
    void OutputPymolLabelProteinModelFromString< t_PrecisionType>::ListNodeMembers
    (
      const Node< std::string, t_PrecisionType> &NODE,
      const std::string &LIST_FILENAME,
      const std::string &PREFIX,
      const std::string &POSTFIX
    )
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, LIST_FILENAME);

      // iterate through the members of NODE in order to list the members
      for
      (
        util::SiPtrList< std::string>::const_iterator
          member_itr( NODE.GetMembers().Begin()), member_itr_end( NODE.GetMembers().End());
        member_itr != member_itr_end;
        ++member_itr
      )
      {
        write << PREFIX + ( **member_itr) + POSTFIX << '\n';
      }

      io::File::CloseClearFStream( write);
    }

    //! @brief removes characters from string which could cause problems within python code
    //! @param STRING the string which will be cleaned
    template< typename t_PrecisionType>
    void OutputPymolLabelProteinModelFromString< t_PrecisionType>::CleanString( std::string &STRING)
    {
      {
        const util::StringReplacement replacer( util::StringReplacement::e_Any, "/", "X");
        replacer.ReplaceEachIn( STRING);
      }
      {
        const util::StringReplacement replacer( util::StringReplacement::e_Any, "-", "X");
        replacer.ReplaceEachIn( STRING);
      }
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_PYMOL_LABEL_PROTEIN_MODEL_FROM_STRING_H_ 
