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

#ifndef BCL_BIOL_SASA_POINT_H_
#define BCL_BIOL_SASA_POINT_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasaPoint
    //! @brief Wrapper Class of SasaData
    //! @details Class to explicitly list Atom number, Solvent Excluded Surface, and Solvent Accessible Surface      //!
    //! @see @link example_biol_SasaPoint.cpp @endlink
    //! @author putnamdk
    //! @date April 9, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasaPoint :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief AtomNumber represents the atom order in the PDB ( zero based indexing)
      size_t m_AtomNumber;

      //! @brief Solvent Excluded Surface represents the topological boundary of the union of all possible probes
      //! @brief which do not overlap with the molecule.  The number is stored as a percentage
      double m_SolventExcludedSurface;

      //! @brief Solvent Accessible Surface is traced out by the center of the probe representing a solvent molecule
      //! @brief The number is stored as a percentage
      double m_SolventAccessibleSurface;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasaPoint();

      //! @brief constructor from given input data
      //! @param ATOMNUM - The atom number in the PDB
      //! @param SES     - The Solvent Excluded Surface
      //! @param SAS     - The Solvent Accessible Surface
      SasaPoint( const size_t ATOMNUM, const double SES, const double SAS);

      //! @brief Clone function
      //! @return pointer to new SasaPoint
      SasaPoint *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Accessor Function to Private data variable
      //! @return the Atom number value
      const size_t &GetAtomNumber() const
      {
        return m_AtomNumber;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The Solvent Excluded Surface (SES)
      const double &GetSolventExcludedSurface() const
      {
        return m_SolventExcludedSurface;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The Solvent Accessible Surface (SAS)
      const double &GetSolventAccessibleSurface() const
      {
        return m_SolventAccessibleSurface;
      }

    ////////////////
    // operators  //
    ////////////////

      //! @brief compare two SasaPoint
      //! @param POINT the point to compare to this point
      //! @return if *this and POINT have identical data
      bool operator ==( const SasaPoint &POINT) const
      {
        return
            m_AtomNumber == POINT.m_AtomNumber &&
            m_SolventExcludedSurface == POINT.m_SolventExcludedSurface &&
            m_SolventAccessibleSurface == POINT.m_SolventAccessibleSurface;
      }

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

    }; // class SasaPoint
  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_SASA_POINT_H_
