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

#ifndef BCL_SDF_MDL_HEADER_H_
#define BCL_SDF_MDL_HEADER_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MdlHeader
    //! @brief manages all MDL lines from a MDL file/block in terms of reading and parsing mdl lines of
    //!        certain line types as well as writing the lines back to create a MDL file from a MDL Handler.
    //!
    //! @see @link example_sdf_mdl_header.cpp @endlink
    //! @author mendenjl
    //! @date Feb 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MdlHeader :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t m_NumberAtoms; //!< Number of atoms declared in the mdl
      size_t m_NumberBonds; //!< Number of bonds declared in the mdl
      bool   m_IsValid;     //!< Whether reading of the mdl should proceed

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MdlHeader();

      //! @brief constructor from # atoms, # bonds
      MdlHeader( const size_t &NUMBER_ATOMS, const size_t &NUMBER_BONDS);

      //! @brief Clone function
      //! @return pointer to new MdlHeader
      MdlHeader *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get # atom lines declared
      //! @return # atom lines declared
      const size_t &GetNumberAtoms() const
      {
        return m_NumberAtoms;
      }

      //! @brief get # bond lines declared
      //! @return # bond lines declared
      const size_t &GetNumberBonds() const
      {
        return m_NumberBonds;
      }

      //! @brief return the string of the header
      //! @return the string of the header
      std::string ToMdlLine() const;

      //! @brief return the string of the header
      //! @param MDL_HEADER the header line
      //! @param NUMBER_DESC_LINES the number of description lines read so far. If there
      //! @return the string of the header
      void SetFromMdlLine( const std::string &MDL_HEADER, const size_t NUMBER_DESC_LINES);

      //! @brief check if the handler is valid
      //! @return true if the mdl handler is in a valid state
      bool IsValid() const
      {
        return m_IsValid;
      }

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MdlHeader

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_MDL_HEADER_H_

