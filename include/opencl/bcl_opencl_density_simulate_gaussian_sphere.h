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

#ifndef BCL_OPENCL_DENSITY_SIMULATE_GAUSSIAN_SPHERE_H_
#define BCL_OPENCL_DENSITY_SIMULATE_GAUSSIAN_SPHERE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_vector.h"
#include "density/bcl_density_simulate_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DensitySimulateGaussianSphere
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_opencl_DensitySimulateGaussianSphere.cpp @endlink
    //! @author woetzen
    //! @date Nov 30, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DensitySimulateGaussianSphere :
      public density::SimulateInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_GridSpacing; //!< grid spacing for the simulated density map
      double         m_Resolution;  //!< resolution for the simulated density map
      size_t         m_Margin;      //!< margin to be added to grid extent as determined by atom coordinates

      CommandQueue   m_CommandQueue; //!< the command queue
      cl::Program    m_Program;      //!< the opencl program

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DensitySimulateGaussianSphere();

      //! @brief Clone function
      //! @return pointer to new DensitySimulateGaussianSphere
      DensitySimulateGaussianSphere *Clone() const;

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

      //! @brief is this class compatible with given command queue
      //! @param COMMAND_QUEUE the command queue this object would operate on
      //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
      bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const;

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE);

      //! @brief simulate a density map for a buffer of atoms into given grid dimensions
      //! @param ATOMS a Buffer of atoms with their weight
      //! @param NR_ATOMS number of atoms in Buffer with size 4 * NR_ATOMS
      //! @param INDEX the index of the grid
      //! @param DIMENSION the dimension of the grid
      Vector< double> Simulate
      (
        const Matrix< double> &ATOMS,
        const size_t NR_ATOMS,
        const storage::VectorND< 3, int> &INDEX,
        const storage::VectorND< 3, size_t> &DIMENSIONS
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief generate simulated density from given list of atoms
      //! @param ATOMS siptrvector of atoms
      //! @return a simulated density map
      density::Map operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const;

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
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief atoms to device buffer
      //! @param ATOMS siptrvector of atoms
      //! @return Buffer a buffer on the device that is associated with this command queue
      Matrix< double> AtomsToDevice( const util::SiPtrVector< const biol::Atom> &ATOMS) const;

      //! @brief round up dimensions so that they are compatible with work group size
      //! @param DIMENSIONS the original dimensions
      //! @return storage::VectorND< 3, size_t> the new diemsions rounded up
      storage::VectorND< 3, size_t> RoundUpDimensions( const storage::VectorND< 3, size_t> &DIMENSIONS) const;

    private:

      //! @brief convert atoms to padded matrix
      //! @param ATOMS siptrvector of atoms
      //! @return linal::Matrix< double> a matrix with 4 cols (3 coordinates with 1 weight) and number of atoms + x rows
      //!         ( which have the coordinates of the last row with 0 weight) so that they are a multiple of 16
      static linal::Matrix< double> AtomsToMatrix( const util::SiPtrVector< const biol::Atom> &ATOMS, const double BLOB_C);

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occurred, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION);

    }; // class DensitySimulateGaussianSphere

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_DENSITY_SIMULATE_GAUSSIAN_SPHERE_H_ 
