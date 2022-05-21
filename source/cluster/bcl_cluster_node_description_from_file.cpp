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

// include forward header of this class
#include "cluster/bcl_cluster_node_description_from_file.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NodeDescriptionFromFile::s_Instance
    (
      GetObjectInstances().AddInstance( new NodeDescriptionFromFile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NodeDescriptionFromFile::NodeDescriptionFromFile() :
      m_Descriptions()
    {
    }

    //! @brief constructor taking a filename and filling "m_Descriptions" based on the file's content
    //! @param FILENAME the path and name of the file which contains the strings and their descriptions
    NodeDescriptionFromFile::NodeDescriptionFromFile( const std::string &FILENAME) :
      m_Descriptions( GetDescriptions( FILENAME))
    {
    }

    //! @brief copy constructor
    //! @param DESCRIPTION the object which will be copied into this instance
    NodeDescriptionFromFile::NodeDescriptionFromFile
    (
      const NodeDescriptionFromFile &DESCRIPTION
    ) :
      m_Descriptions( DESCRIPTION.m_Descriptions)
    {
    }

    //! @brief virtual copy constructor
    //! @return new copy of this class
    NodeDescriptionFromFile *NodeDescriptionFromFile::Clone() const
    {
      return new NodeDescriptionFromFile( *this);
    }

    //! @brief virtual destructor
    NodeDescriptionFromFile::~NodeDescriptionFromFile()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &NodeDescriptionFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking a Wrapper<string> and returns a double which is the desription of the object
    //! @param MEMBER the object which was clustered for which a description is desired
    //! @return returns a double which is the desription of "MEMBER"
    double NodeDescriptionFromFile::operator()( const std::string &MEMBER) const
    {
      // create const iterator "itr" to the "MEMBER" key in "m_Descriptions"
      storage::Map< std::string, double>::const_iterator itr( m_Descriptions.Find( MEMBER));

      // make sure that "itr" is not at the end of "m_Descriptions" - "MEMBER" could not be found in "m_Descriptions"
      if( itr == m_Descriptions.End())
      {
        // write a warning message that "MEMBER" could not be found
        BCL_MessageStd( MEMBER + " could not be found in descriptor list");

        // return undefined double
        return util::GetUndefined< double>();
      }

      // return the value description that corresponds to "MEMBER"
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NodeDescriptionFromFile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Descriptions, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &NodeDescriptionFromFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Descriptions, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief GetDescriptions takes a file and reads in the descriptions for the clustered objects
    //! @param "FILENAME" the path and name of the file which contains the strings and their descriptions
    //! @return return a map which contains the clustered string objects and each numberical description
    storage::Map< std::string, double> NodeDescriptionFromFile::GetDescriptions
    (
      const std::string &FILENAME
    ) const
    {
      // create IFStream "read"
      io::IFStream read;

      // open "read" and bind it to "FILENAME"
      io::File::MustOpenIFStream( read, FILENAME);

      // create Vector of Vectors of strings "split_lines" and initialize with the split lines in "FILENAME"
      const storage::Vector< storage::Vector< std::string> > split_lines
      (
        util::SplittedStringLineListFromIStream( read)
      );

      // create Map "descriptions" which will store all of the descriptions and clustered strings in "FILENAME"
      storage::Map< std::string, double> descriptions;

      // iterate through "split_lines" in order to fill "descriptions"
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr( split_lines.Begin()), itr_end( split_lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // create const size_t "line_num_cols" and initialize with the number of columns in
        // the current line denoted by "itr"
        const size_t line_num_cols( itr->GetSize());

        // true if there are less than two columns in the current line of "FILENAME" denoted by "itr"
        if( line_num_cols < 2)
        {
          // write a warning message that there is not at least 2 columns in current line
          BCL_MessageStd
          (
            "Only " + util::Format()( line_num_cols) + " on line " +
            util::Format()( itr - split_lines.Begin()) + " of file " + FILENAME
          );

          // go to next line
          continue;
        }

        // true if there are more than two columns in the current line of "FILENAME" denoted by "itr"
        // any addition columns will be ignored
        else if( line_num_cols > 2)
        {
          // write a warning message that additional columns will be ignored
          BCL_MessageStd
          (
            "More than 2 columns. Additional columns in " + FILENAME + " will be ignored"
          );
        };

        // create const reference to string "string_object" and initialize with the string object denoted by "itr"
        const std::string &string_object( itr->operator()( 0));

        // create const double "string_object" and initialize with the double description denoted by "itr"
        const double double_description( util::ConvertStringToNumericalValue< double>( itr->operator()( 1)));

        // create pair "insertion" and initialize with the results of trying to
        // insert the string and its numerical description into "descriptions" from the first and second columns
        std::pair< storage::Map< std::string, double>::const_iterator, bool> insertion
        (
          descriptions.Insert
          (
            std::pair< std::string, double>
            (
              string_object, double_description
            )
          )
        );

        // true if the insertion of the string and its description was unsuccessful
        if( !insertion.second)
        {
          // write warning message that the insertion failed maybe because that key already exists in "descriptions"
          BCL_MessageStd
          (
            "Failed to insert " + string_object + " and its description " +
            util::Format()( double_description) + " Maybe it already has a description"
          );
        }
      }

      // return "descriptions"
      return descriptions;
    }

  } // namespace cluster
} // namespace bcl
