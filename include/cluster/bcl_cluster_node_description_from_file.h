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

#ifndef BCL_CLUSTER_NODE_DESCRIPTION_FROM_FILE_H_
#define BCL_CLUSTER_NODE_DESCRIPTION_FROM_FILE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NodeDescriptionFromFile
    //! @brief NodeDescriptionFromFile is for getting the numerical descriptions of clustered objects
    //! @details This class is specifically for when strings are clustered such as when clustering from a matrix of RMSDs and
    //! corresponding file names
    //! The describing number could be anything about the object, for example, its score or its RMSD to a native
    //! structure
    //! It reads in a file and uses the first two columns to associate a string object with its numerical description
    //! If an object occurs more than once, only the first description will be kept
    //! If more than 2 columns occur in a line in the file they will be ignored
    //!
    //! @see @link example_cluster_node_description_from_file.cpp @endlink
    //! @author alexanns
    //! @date September 9, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NodeDescriptionFromFile :
      public util::FunctionInterface< std::string, double>
    {
    protected:

    //////////
    // data //
    //////////

      //! "m_Descriptions" stores the clustered objects and their corresponding description
      storage::Map< std::string, double> m_Descriptions;

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
      NodeDescriptionFromFile();

      //! @brief constructor taking a filename and filling "m_Descriptions" based on the file's content
      //! @param FILENAME the path and name of the file which contains the strings and their descriptions
      explicit NodeDescriptionFromFile( const std::string &FILENAME);

      //! @brief copy constructor
      //! @param DESCRIPTION the object which will be copied into this instance
      NodeDescriptionFromFile( const NodeDescriptionFromFile &DESCRIPTION);

      //! @brief virtual copy constructor
      //! @return new copy of this class
      NodeDescriptionFromFile *Clone() const;

      //! @brief virtual destructor
      virtual ~NodeDescriptionFromFile();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking a Wrapper<string> and returns a double which is the desription of the object
      //! @param MEMBER the object which was clustered for which a description is desired
      //! @return returns a double which is the desription of "MEMBER"
      double operator()( const std::string &MEMBER) const;

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

    protected:

      //! @brief GetDescriptions takes a file and reads in the descriptions for the clustered objects
      //! @param "FILENAME" the path and name of the file which contains the strings and their descriptions
      //! @return return a map which contains the clustered string objects and each numberical description
      storage::Map< std::string, double> GetDescriptions( const std::string &FILENAME) const;

    }; // class NodeDescriptionFromFile

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_NODE_DESCRIPTION_FROM_FILE_H_ 
