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

#ifndef BCL_SDF_BOND_INFO_H_
#define BCL_SDF_BOND_INFO_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_types.h"
#include "chemistry/bcl_chemistry_constitutional_bond_types.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondInfo
    //! @brief Manages generic bond info; some of this bond info may be converted to/from bond lines from the mdl,
    //!        however, the data structure and members are independent of the mdl format
    //!
    //! @see @link example_sdf_bond_info.cpp @endlink
    //! @author mendenjl
    //! @date Mar 02, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondInfo :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t                             m_AtomIndexLow;   //!< Smaller atom index
      size_t                             m_AtomIndexHigh;  //!< Larger atom index
      chemistry::ConfigurationalBondType m_BondType;       //!< Bond type

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BondInfo();

      //! @brief constructor from members
      BondInfo
      (
        const size_t &ATOM_INDEX_A,
        const size_t &ATOM_INDEX_B,
        const chemistry::ConfigurationalBondType &BOND_TYPE
      );

      //! @brief constructor from members
      BondInfo
      (
        const size_t &ATOM_INDEX_A,
        const size_t &ATOM_INDEX_B,
        const chemistry::ConstitutionalBondType &BOND_TYPE
      );

      //! @brief Clone function
      //! @return pointer to new BondInfo
      BondInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the smaller of the two atom indices
      //! @return the smaller of the two atom indices
      const size_t &GetAtomIndexLow() const
      {
        return m_AtomIndexLow;
      }

      //! @brief get the larger of the two atom indices
      //! @return the larger of the two atom indices
      const size_t &GetAtomIndexHigh() const
      {
        return m_AtomIndexHigh;
      }

      //! @brief get the constitutional bond type
      //! @return the constitutional bond type
      const chemistry::ConstitutionalBondType &GetConstitutionalBondType() const
      {
        return m_BondType->GetConstitutionalBondType();
      }

      //! @brief get the configurational bond type
      //! @return the configurational bond type
      const chemistry::ConfigurationalBondType &GetConfigurationalBondType() const
      {
        return m_BondType;
      }

      //! @brief check if the handler is valid
      //! @return true if the mdl handler is in a valid state
      bool IsValid() const
      {
        return m_BondType.IsDefined();
      }

      //! @brief set the chirality from a given property string
      //! @param ISOMETRY the desired isometry
      void SetIsometry( const chemistry::BondIsometry &ISOMETRY);

    ////////////////
    // operations //
    ////////////////

      //! @brief determine if the specified line fits the formatting for a bond line
      //! @param LINE the line to check
      //! @return true if the line fits the formatting
      static bool FormattedAsMdlBondLine( const std::string &LINE);

      //! @brief extract information from a line of the sdf believed to be an mdl bond line
      //! @param LINE line from an sdf file that is believed to be a bond line
      //! @return reference to this
      //! This only extracts the bond order or whether the bond is aromatic
      //! the bond type can be refined by later calls to set functions
      BondInfo &ExtractMdlBondLineInfo( const std::string &LINE);

      //! @brief create an mdl bond line string out of the bond info
      //! @return an mdl bond line string created from the bond info
      std::string ToMdlBondLine() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief compare this bond info to other bond info
      //! @param BOND_INFO the object to compare this bond info to
      //! @return true if the bond info is the same
      bool operator ==( const BondInfo &BOND_INFO) const;

      //! @brief compare this bond info to other bond info
      //! @param BOND_INFO the object to compare this bond info to
      //! @return true if the bond info is less than BOND_INFO
      bool operator <( const BondInfo &BOND_INFO) const;

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

    }; // class BondInfo

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_BOND_INFO_H_

