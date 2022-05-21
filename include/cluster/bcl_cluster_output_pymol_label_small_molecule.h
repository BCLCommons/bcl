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

#ifndef BCL_CLUSTER_OUTPUT_PYMOL_LABEL_SMALL_MOLECULE_H_
#define BCL_CLUSTER_OUTPUT_PYMOL_LABEL_SMALL_MOLECULE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_pymol_label_string.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "math/bcl_math.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputPymolLabelSmallMolecule
    //! @brief For use with the OutputPymol object and is an object that can
    //! create labels for nodes when the object being clustering is a string that denotes an sdf file for some
    //! small molecules
    //! @details It is derived from util::FunctionInterface and takes a storage::Triplet and returns a std::string which
    //! is the text necessary to put into a Python script that can be run by Pymol in order to generate the label
    //! The label consists of the
    //! actual molecule which is the node center, the number of members in that node, and the girth of the node
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
    //! Molecules must be numbered starting with zero in the table ( 0 to N - 1, with N begin the number of molecules)
    //!
    //! @see @link example_cluster_output_pymol_label_small_molecule.cpp @endlink
    //! @author alexanns
    //! @date Dec 2, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_PrecisionType>
    class OutputPymolLabelSmallMolecule :
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

    private:

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

      //! holds all of the small molecules from an sdf file
      util::ShPtrVector< sdf::MdlHandler> m_MdlHandlers;

      //! color for lables
      util::Color m_Color;

      //! true if each member of every labeled node should be written to an mdl file
      bool m_WriteMembersMdlFiles;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      OutputPymolLabelSmallMolecule();

      //! @brief constructor taking parameters for each member variable
      //! @param MEMBER_DISTANCE_FUNCTION the function used to compare the members and get the distance between members
      //!        within the dendrogram
      //! @param DISTANCE_COMPARATOR the way in which the distance calculated by "m_MemberDistanceFunction" should
      //!        be compared. For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
      //! @param SDF_FILENAME the filename of the sdf file containing the small molecules
      //! @param COLOR the color for the labels
      //! @param WRITE_MEMBER_MDL true if each member of every labeled node should be written to an mdl file
      OutputPymolLabelSmallMolecule
      (
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
        > &MEMBER_DISTANCE_FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &DISTANCE_COMPARATOR,
        const std::string &SDF_FILENAME,
        const util::Color &COLOR,
        const bool WRITE_MEMBER_MDL = false
      );

      //! @brief Clone function
      //! @return pointer to new OutputPymolSmallMolecule
      OutputPymolLabelSmallMolecule *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief provides the text necessary for loading an object into pymol with a given pymol object name
      //! @param FILENAME the filename where the object will be loaded from
      //! @param PYMOL_OBJECT_NAME the name of the object as it will be loaded into pymol as
      //! @return the string which contains the text that can go into a python script to load the file into pymol
      static std::string LoadPymolObject
      (
        const std::string &FILENAME, const std::string &PYMOL_OBJECT_NAME
      );

      //! @brief provides the text necessary for translating an object in pymol to a new location
      //! @param X_COORDINATE the desired x translation of the object
      //! @param Y_COORDINATE the desired y translation of the object
      //! @param Z_COORDINATE the desired z translation of the object
      //! @param PYMOL_OBJ_NAME the name of the object in pymol that will be translated
      //! @return the string which has the text necessary for translating the desired object the desired amount
      static std::string TranslatePymolObject
      (
        const double X_COORDINATE, const double Y_COORDINATE, const double Z_COORDINATE,
        const std::string &PYMOL_OBJ_NAME
      );

      //! @brief gives the text necessary to create a label based on the object of a cluster
      //! @param X_COORDINATE the x coordinate where the label will be placed
      //! @param Y_COORDINATE_TOP the top coordinate of the node cylinder which the label is for
      //! @param Y_COORDINATE_BOTTOM the bottom coordinate of the node cylinder which the label is for
      //! @param NODE_OBJECT the object which is the center of the node which is going to be the label
      //! @param NODE_CYLINDER_DIAMETER the diameter of the node which is going to be labeled
      //! @pararm COLOR color for label
      //! @return a string which has all the text necessary to create a label for a node about the filename
      static std::string GetObjectNameLabel
      (
        const double X_COORDINATE, const double Y_COORDINATE_TOP,
        const double Y_COORDINATE_BOTTOM, const util::SiPtr< const std::string> NODE_CENTER,
        const double NODE_CYLINDER_DIAMETER,
        const util::Color &COLOR
      );

    protected:

      //! @brief writes all the members of the given node to individual mdl files
      //! @param NODE the node whose members will be written into mdl files
      //! @param PATH the path where the mdl files will be written to
      void WriteNodeMemberMdlFiles( const Node< std::string, t_PrecisionType> &NODE, const std::string &PATH) const;

      //! @brief writes an mdl file for a give node member
      //! @param FILENAME the name of the mdl file that will be written to
      //! @param MEMBER the member which will be written to the mdl file
      //! @return the Mdl handler which is used to write out the mdl file
      const sdf::MdlHandler &WriteMDLFile( const std::string &FILENAME, const std::string &MEMBER) const;

    }; // class OutputPymolLabelSmallMolecule

    //! single instance of that class
    template< typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> OutputPymolLabelSmallMolecule< t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new OutputPymolLabelSmallMolecule< t_PrecisionType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_PrecisionType>
    OutputPymolLabelSmallMolecule< t_PrecisionType>::OutputPymolLabelSmallMolecule() :
      m_MemberDistanceFunction(),
      m_DistanceComparator(),
      m_MdlHandlers(),
      m_Color( util::GetColors().e_Magenta),
      m_WriteMembersMdlFiles( false)
    {
    }

    //! @brief constructor taking parameters for each member variable
    //! @param MEMBER_DISTANCE_FUNCTION the function used to compare the members and get the distance between members
    //!        within the dendrogram
    //! @param DISTANCE_COMPARATOR the way in which the distance calculated by "m_MemberDistanceFunction" should
    //!        be compared. For example, "<" for RMSD (a dissimilarity measure) or ">" for a similarity measure
    //! @param SDF_FILENAME the filename of the sdf file containing the small molecules
    //! @param COLOR the color for the labels
    //! @param WRITE_MEMBER_MDL true if each member of every labeled node should be written to an mdl file
    template< typename t_PrecisionType>
    OutputPymolLabelSmallMolecule< t_PrecisionType>::OutputPymolLabelSmallMolecule
    (
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const std::string> >, t_PrecisionType>
      > &MEMBER_DISTANCE_FUNCTION,
      const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &DISTANCE_COMPARATOR,
      const std::string &SDF_FILENAME,
      const util::Color &COLOR,
      const bool WRITE_MEMBER_MDL
    ) :
      m_MemberDistanceFunction( MEMBER_DISTANCE_FUNCTION),
      m_DistanceComparator( DISTANCE_COMPARATOR),
      m_MdlHandlers(),
      m_Color( COLOR),
      m_WriteMembersMdlFiles( WRITE_MEMBER_MDL)
    {
      io::IFStream read;
      io::File::MustOpenIFStream( read, SDF_FILENAME);

      do
      {
        // create new handler, which should read to the next terminator
        util::ShPtr< sdf::MdlHandler> sp_current_handler( new sdf::MdlHandler( read));

        // check if the mdl was successfully read
        if( sp_current_handler->IsValid())
        {
          m_MdlHandlers.PushBack( sp_current_handler);
        }
        else
        {
          break;
        }
      } while( read.good() && !read.eof());
    }

    //! @brief Clone function
    //! @return pointer to new OutputPymolLabelSmallMolecule
    template< typename t_PrecisionType>
    OutputPymolLabelSmallMolecule< t_PrecisionType> *OutputPymolLabelSmallMolecule< t_PrecisionType>::Clone() const
    {
      return new OutputPymolLabelSmallMolecule< t_PrecisionType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_PrecisionType>
    const std::string &OutputPymolLabelSmallMolecule< t_PrecisionType>::GetClassIdentifier() const
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
    std::string OutputPymolLabelSmallMolecule< t_PrecisionType>::operator()
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

      BCL_MessageDbg( "center is " + util::Format()( node_center));

      // get the name of the mol file that will be written out for the node center
      const std::string mol_filename( io::File::SplitToPathAndFileName( script_filename).First() + ( *node_center) + ".mol");

      // write the center to an mdl file
      const sdf::MdlHandler &handler( WriteMDLFile( mol_filename, *node_center));

      std::string label_command;

      std::string mol_pymol_obj_name( ( *node_center) + "_" + util::Format()( node->GetIdentifier()));
      util::StringReplacement replacer( util::StringReplacement::e_Any, "/", "_");
      replacer.ReplaceEachIn( mol_pymol_obj_name);

      // load the mdl file into pymol
      label_command += LoadPymolObject( mol_filename, mol_pymol_obj_name);

      // create the small molecule from the handler
      chemistry::FragmentComplete small_mol( sdf::FragmentFactory::MakeFragment( handler));

      // size of molecule
      const double volume( descriptor::GetCheminfoProperties().calc_VdwVolume->SumOverObject( small_mol)( 0));
      BCL_MessageDbg( "the volume of the molecule is " + util::Format()( volume));
      const double radius( math::Pow( ( 3.0 / 4.0) * ( 1.0 / math::g_Pi), 1.0 / 3.0));
      BCL_MessageDbg( "the radius of the molecule is " + util::Format()( radius));

      // move the molecule to the position of the node
      const double translate_x_coordinate( x_coordinate);
      const double translate_y_coordinate( ( y_coordinate_top + y_coordinate_bottom) / 2.0);
      const double translate_z_coordinate( 20.0 * radius + node_cylinder_diameter / 2.0);
      label_command += TranslatePymolObject
        ( translate_x_coordinate, translate_y_coordinate, translate_z_coordinate, mol_pymol_obj_name);

      // add command to show the size and girth of "node"
      label_command += OutputPymolLabelString< t_PrecisionType>::template GetLabelCommand< std::string>
      (
        *node, "size_girth_" + mol_pymol_obj_name,
        x_coordinate, y_coordinate_top, y_coordinate_bottom, node_cylinder_diameter, m_Color
      );

      // get a label to label the filename of the molecule which is the center of the cluster
      label_command += GetObjectNameLabel
        ( x_coordinate, y_coordinate_top, y_coordinate_bottom, node_center, node_cylinder_diameter, m_Color);

      // true if write the members of the node
      if( m_WriteMembersMdlFiles)
      {
        WriteNodeMemberMdlFiles( *node, io::File::SplitToPathAndFileName( script_filename).First());
      }

      return label_command;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_PrecisionType>
    std::istream &OutputPymolLabelSmallMolecule< t_PrecisionType>::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_MemberDistanceFunction, ISTREAM);
      io::Serialize::Read( m_DistanceComparator    , ISTREAM);
      io::Serialize::Read( m_Color                 , ISTREAM);
      io::Serialize::Read( m_WriteMembersMdlFiles  , ISTREAM);

      do
      {
        // create new handler, which should read to the next terminator
        util::ShPtr< sdf::MdlHandler> sp_current_handler( new sdf::MdlHandler( ISTREAM));

        // check if the mdl was successfully read
        if( sp_current_handler->IsValid())
        {
          m_MdlHandlers.PushBack( sp_current_handler);
        }
        else
        {
          break;
        }
      } while( ISTREAM.good() && !ISTREAM.eof());

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_PrecisionType>
    std::ostream &OutputPymolLabelSmallMolecule< t_PrecisionType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      // write members
      io::Serialize::Write( m_MemberDistanceFunction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceComparator    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Color                 , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WriteMembersMdlFiles  , OSTREAM, INDENT) << '\n';

      // write the mdl lines
      for
      (
        util::ShPtrVector< sdf::MdlHandler>::const_iterator
          mdl_itr( m_MdlHandlers.Begin()), mdl_itr_end( m_MdlHandlers.End());
        mdl_itr != mdl_itr_end;
        ++mdl_itr
      )
      {
        ( *mdl_itr)->WriteToSDF( OSTREAM);
      }

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief provides the text necessary for loading an object into pymol with a given pymol object name
    //! @param FILENAME the filename where the object will be loaded from
    //! @param PYMOL_OBJECT_NAME the name of the object as it will be loaded into pymol as
    //! @return the string which contains the text that can go into a python script to load the file into pymol
    template< typename t_PrecisionType>
    std::string OutputPymolLabelSmallMolecule< t_PrecisionType>::LoadPymolObject
    (
      const std::string &FILENAME, const std::string &PYMOL_OBJECT_NAME
    )
    {
      return "cmd.load( \"" + FILENAME + "\", \"" + PYMOL_OBJECT_NAME + "\")\n";
    }

    //! @brief provides the text necessary for translating an object in pymol to a new location
    //! @param X_COORDINATE the desired x translation of the object
    //! @param Y_COORDINATE the desired y translation of the object
    //! @param Z_COORDINATE the desired z translation of the object
    //! @param PYMOL_OBJ_NAME the name of the object in pymol that will be translated
    //! @return the string which has the text necessary for translating the desired object the desired amount
    template< typename t_PrecisionType>
    std::string OutputPymolLabelSmallMolecule< t_PrecisionType>::TranslatePymolObject
    (
      const double X_COORDINATE, const double Y_COORDINATE, const double Z_COORDINATE, const std::string &PYMOL_OBJ_NAME
    )
    {
      return "cmd.translate( [" + util::Format()( X_COORDINATE) + ", " + util::Format()( Y_COORDINATE) + ", " + util::Format()( Z_COORDINATE) + "], \"" + PYMOL_OBJ_NAME + "\")\n";
    }

    //! @brief gives the text necessary to create a label based on the object of a cluster
    //! @param X_COORDINATE the x coordinate where the label will be placed
    //! @param Y_COORDINATE_TOP the top coordinate of the node cylinder which the label is for
    //! @param Y_COORDINATE_BOTTOM the bottom coordinate of the node cylinder which the label is for
    //! @param NODE_OBJECT the object which is the center of the node which is going to be the label
    //! @param NODE_CYLINDER_DIAMETER the diameter of the node which is going to be labeled
    //! @pararm COLOR color for label
    //! @return a string which has all the text necessary to create a label for a node about the filename
    template< typename t_PrecisionType>
    std::string OutputPymolLabelSmallMolecule< t_PrecisionType>::GetObjectNameLabel
    (
      const double X_COORDINATE, const double Y_COORDINATE_TOP,
      const double Y_COORDINATE_BOTTOM, const util::SiPtr< const std::string> NODE_CENTER,
      const double NODE_CYLINDER_DIAMETER,
      const util::Color &COLOR
    )
    {
      // name of variable in python code
      std::string object_name( "center_" + io::File::RemoveFullExtension( *NODE_CENTER));

      OutputPymolLabelProteinModelFromString< t_PrecisionType>::CleanString( object_name);

      // add comand to create Python object "file_object_name"
      std::string label( object_name + "=[]\n");

      // create const double "text_block_size" and initialize with the results of GetTextBlockSize
      const double text_block_size
      (
        OutputPymolLabelString< t_PrecisionType>::GetTextBlockSize( Y_COORDINATE_TOP, Y_COORDINATE_BOTTOM)
      );

      // create const double "text_cylinder_radius" and initialize with the radius that the cylinders making up
      // the text should have
      const double text_cylinder_radius( OutputPymolLabelString< t_PrecisionType>::GetCylinderTextRadius( text_block_size));

      // add command to create Python code to create Pymol object which is the text of "node" center string which is
      // a pdb file
      label += OutputPymolLabelString< t_PrecisionType>::GetCylinderText
      (
        object_name,
        X_COORDINATE,
        ( ( ( Y_COORDINATE_TOP + Y_COORDINATE_BOTTOM) / 2.0) + Y_COORDINATE_TOP) / 2.0,
        OutputPymolLabelString< t_PrecisionType>::GetTextZCoordinate( NODE_CYLINDER_DIAMETER, text_cylinder_radius),
        *NODE_CENTER, text_cylinder_radius, text_block_size, COLOR
      );

      // add command to load the label "file_object_name" object in Pymol
      label += "labels_all.extend( " + object_name + ")\n";

      return label;
    }

    //! @brief writes all the members of the given node to individual mdl files
    //! @param NODE the node whose members will be written into mdl files
    //! @param PATH the path where the mdl files will be written to
    template< typename t_PrecisionType>
    void OutputPymolLabelSmallMolecule< t_PrecisionType>::WriteNodeMemberMdlFiles
    (
      const Node< std::string, t_PrecisionType> &NODE,
      const std::string &PATH
    ) const
    {
      // iterate through the members of NODE in order to write out mdl files for them
      for
      (
        util::SiPtrList< std::string>::const_iterator
          member_itr( NODE.GetMembers().Begin()), member_itr_end( NODE.GetMembers().End());
        member_itr != member_itr_end;
        ++member_itr
      )
      {
        // construct filename for mdl
        std::string mdl_filename( PATH + "node_" + util::Format()( NODE.GetIdentifier()) + "_" + **member_itr + ".mol");

        // write the mdl for the current member to the file
        WriteMDLFile( mdl_filename, **member_itr);
      }
    }

    //! @brief writes an mdl file for a give node member
    //! @param FILENAME the name of the mdl file that will be written to
    //! @param MEMBER the member which will be written to the mdl file
    //! @return the Mdl handler which is used to write out the mdl file
    template< typename t_PrecisionType>
    const sdf::MdlHandler &OutputPymolLabelSmallMolecule< t_PrecisionType>::WriteMDLFile
    (
      const std::string &FILENAME, const std::string &MEMBER
    ) const
    {
      BCL_MessageDbg( "Writing mdl file " + FILENAME);

      size_t index( util::ConvertStringToNumericalValue< size_t>( MEMBER));
      BCL_Assert
      (
        index <= m_MdlHandlers.GetSize(),
        "Given molecule file for output_pymol_label_small_molecule had fewer molecules (" +
        util::Format()( m_MdlHandlers.GetSize()) + " than the distance matrix (which was of size at least "
        + util::Format()( index) + ")"
      );

      // write the mdl file
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);
      m_MdlHandlers( index)->WriteToSDF( write);
      io::File::CloseClearFStream( write);

      return *m_MdlHandlers( index);
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_PYMOL_LABEL_SMALL_MOLECULE_H_
