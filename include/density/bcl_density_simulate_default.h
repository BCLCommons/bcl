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

#ifndef BCL_DENSITY_SIMULATE_DEFAULT_H_
#define BCL_DENSITY_SIMULATE_DEFAULT_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_density_simulate_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SimulateDefault
    //! @brief this class simulates a density map from a given set of atoms
    //!
    //! @see @link example_density_simulate_default.cpp @endlink
    //! @author woetzen
    //! @date Jun 21, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SimulateDefault :
      public SimulateInterface
    {

    public:

    ///////////
    // types //
    ///////////

      //! @enum Kernel
      //! @brief can be applied to smooth distributed atom density
      enum Kernel
      {
        e_NoSmoothing,
        e_Gaussian,
        e_Triangular,
        e_SemiEpanechnikov,
        e_Epanechnikov,
        e_HardSphere,
        s_MaxKernel
      };

      //! @brief Kernel as string
      //! @param KERNEL the kernel
      //! @return the string for the kernel
      static const std::string &GetKernelDescriptor( const Kernel &KERNEL)
      {
        static const std::string s_descriptors[] =
        {
          "NoSmoothing",
          "Gaussian",
          "Triangular",
          "SemiEpanechnikov",
          "Epanechnikov",
          "HardSphere",
          GetStaticClassName< Kernel>()
        };
        return s_descriptors[ size_t( KERNEL)];
      }

      //! @brief KernelEnum is used for I/O of smoothing kernel
      typedef util::WrapperEnum< Kernel, &GetKernelDescriptor, s_MaxKernel> KernelEnum;

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_GridSpacing; //!< grid spacing for the simulated density map
      double          m_Resolution;  //!< resolution for the simulated density map
      size_t          m_Margin;      //!< margin to be added to grid extent as determined by atom coordinates
      KernelEnum      m_SmoothingKernel; //!< smoothing kernel to apply

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
      //! @param SMOOTHING_KERNEL kernel for the smoothing
      SimulateDefault( const linal::Vector3D &GRID_SPACING, const double RESOLUTION, const Kernel SMOOTHING_KERNEL);

      //! @brief Clone function
      //! @return pointer to new SimulateDefault
      SimulateDefault *Clone() const;

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

    ///////////////
    // operators //
    ///////////////

      //! @brief generate simulated density from given list of atoms
      //! @param ATOMS siptrvector of atoms
      //! @return a simulated density map
      Map operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const;

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

    }; // class SimulateDefault

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_SIMULATE_DEFAULT_H_
