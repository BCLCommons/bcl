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

#ifndef BCL_DENSITY_SIMULATE_GAUSSIAN_SPHERE_H_
#define BCL_DENSITY_SIMULATE_GAUSSIAN_SPHERE_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_density_simulate_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SimulateGaussianSphere
    //! @brief Simulate a density map from a set of atoms, by representing each atom as a gaussian sphere
    //! @details each atom is represented as a gaussian sphere with a reasonable distance cutoff with the atomic weight
    //!          as intensity and the radius resolution dependent
    //!
    //! @see @link example_density_simulate_gaussian_sphere.cpp @endlink
    //! @author woetzen
    //! @date Jun 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SimulateGaussianSphere :
      public SimulateInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_GridSpacing; //!< grid spacing for the simulated density map
      double         m_Resolution;  //!< resolution for the simulated density map
      size_t         m_Margin;      //!< margin to be added to grid extent as determined by atom coordinates

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from grid spacing and resolution
      //! @param GRID_SPACING the spacing for the density grid
      //! @param RESOLUTION the resolution to simulate for
      SimulateGaussianSphere( const linal::Vector3D &GRID_SPACING, const double RESOLUTION);

      //! @brief Clone function
      //! @return pointer to new SimulateGaussianSphere
      SimulateGaussianSphere *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set the resolution
      //! @param RESOLUTION the resolution for the density map to be generated
      void SetResolution( const double RESOLUTION);

      //! @brief set the resolution
      double GetResolution() const;

      //! @brief set the grid spacing
      //! @param GRID_SPACING the width of a grid element in x, y and z
      void SetGridSpacing( const linal::Vector3D &GRID_SPACING);

      //! @brief set the margin
      //! @param MARGIN number of additional cells next to last atom occupied cells
      void SetMargin( const size_t MARGIN);

    ////////////////
    // operations //
    ////////////////

      //! @brief generate simulated density from given list of atoms
      //! @param ATOMS siptrvector of atoms
      //! @return a simulated density map
      Map operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const;

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SimulateGaussianSphere

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_SIMULATE_GAUSSIAN_SPHERE_H_ 
