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

#ifndef BCL_SDF_ATOM_INFO_H_
#define BCL_SDF_ATOM_INFO_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_types.h"
#include "chemistry/bcl_chemistry_chirality.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomInfo
    //! @brief Manages generic atom info; this atom info may be partially extracted and reconstructed from an MDL line,
    //!        however, the data structure and members are independent of the mdl format
    //!
    //! @see @link example_sdf_atom_info.cpp @endlink
    //! @author mendenjl
    //! @date Mar 02, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomInfo :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      chemistry::AtomType       m_AtomType;      //!< Atom type
      chemistry::ChiralityEnum  m_Chirality;     //!< Chirality
      linal::Vector3D           m_Coordinates;   //!< Coordinates of the atom
      bool                      m_CanAddH;       //!< Whether hydrogens can be added

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AtomInfo();

      //! @brief constructor from all known atom data
      AtomInfo
      (
        const chemistry::AtomType  &ATOM_TYPE,
        const chemistry::Chirality &CHIRALITY,
        const linal::Vector3D      &COORDINATES = linal::Vector3D(),
        const bool                 &CAN_ADD_H = true
      );

      //! @brief Clone function
      //! @return pointer to new AtomInfo
      AtomInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get declared coordinates
      //! @return coordinates
      const linal::Vector3D &GetCoordinates() const
      {
        return m_Coordinates;
      }

      //! @brief get atom type
      //! @return the atom type declared in the sdf file
      const chemistry::AtomType &GetAtomType() const
      {
        return m_AtomType;
      }

      //! @brief get the chirality declared in the sdf file
      //! @return chirality declared in the sdf file
      const chemistry::Chirality &GetChirality() const
      {
        return m_Chirality;
      }

      //! @brief get whether H can be added
      //! @return true if h can be added to this atom
      const bool &CanAddH() const
      {
        return m_CanAddH;
      }

      //! @brief change the charge on the atom type
      //! @param NEW_CHARGE the new charge to give the atom type
      //! This function is necessary if the mdl file contained M  CHG lines that reference this atom
      void SetCharge( const short &CHARGE);

      //! @brief change the chirality of the atom type
      //! @param CHIRALITY the new chirality to give the atom
      //! This function is necessary if the mdl file contained a property called Chirality, one entry for each atom with
      //! 4 bonds
      void SetChirality( const chemistry::Chirality &CHIRALITY);

      //! @brief change the atom type
      //! @param ATOM_TYPE the new atom type
      //! This function is necessary if the mdl file contained a property called AtomType, with one entry for each non-H
      void SetAtomType( const chemistry::AtomType &ATOM_TYPE);

      //! @brief change the coordinates
      //! @param COORDINATES the new coordinates
      void SetCoordinates( const linal::Vector3D &COORDINATES);

    ///////////////
    // operators //
    ///////////////

      //! @brief compare this atom info to other atom info
      //! @param ATOM_INFO the object to compare this atom info to
      //! @return true if the atom info is the same
      bool operator ==( const AtomInfo &ATOM_INFO) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief extract information from the mdl atom line string
      //! @param LINE line from an sdf file that is believed to be an atom line
      //! @return reference to this
      //! This does not extract chirality or the gasteiger atom type (only element type and possibly charge)
      static bool FormattedAsMdlAtomLine( const std::string &LINE);

      //! @brief extract information from the mdl atom line string
      //! @param LINE line from an sdf file that is believed to be an atom line
      //! @return reference to this
      //! This does not extract chirality or the gasteiger atom type (only element type and possibly charge)
      AtomInfo &ExtractMdlAtomLineInfo( const std::string &LINE);

      //! @brief create an mdl atom line string out of the atom info
      //! @return an mdl atom line string created from the atom info
      std::string ToMdlAtomLine() const;

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

    }; // class AtomInfo

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_ATOM_INFO_H_

