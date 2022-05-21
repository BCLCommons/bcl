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

#ifndef BCL_APP_JUFO_H_
#define BCL_APP_JUFO_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Jufo
    //! @brief This is a program that implements the secondary structure prediction method developed by Jens Meiler.
    //! It requires a sequence in fasta file format and a corresponding blast profile in ascii file format and it
    //! returns a jufo output file which list for each residue the likelihood to be in a loop, helix or coil as well
    //! as the one letter code for the most likely option when neighbors are also considered
    //!
    //! @see @link example_app_jufo.cpp @endlink
    //! @author karakam, woetzen, weinerbe
    //! @date Nov 26, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Jufo :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! input fasta file
      util::ShPtr< command::ParameterInterface> m_FastaFileParam;

      //! input ascii file
      util::ShPtr< command::FlagInterface> m_AsciiFileFlag;

      //! output prefix
      util::ShPtr< command::FlagInterface> m_OutputFlag;

      //! use old jufo 3d
      util::ShPtr< command::FlagInterface> m_Jufo3dFlag;

      //! multimer flag
      util::ShPtr< command::FlagInterface> m_MultimerFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      Jufo();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      Jufo *Clone() const
      {
        return new Jufo( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      // instantiate enumerator for Jufo class
      static const ApplicationType Jufo_Instance;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! Main
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      const std::string &GetWebText() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class Jufo

  } // namespace app
} // namespace bcl

#endif // BCL_APP_JUFO_H_
