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

// include header of this class
#include "density/bcl_density_simulate_gaussian_sphere.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SimulateGaussianSphere::s_Instance
    (
      GetObjectInstances().AddInstance( new SimulateGaussianSphere( linal::Vector3D( 2.3, 2.3, 2.3), 6.9))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from grid spacing and resolution
    //! @param GRID_SPACING the spacing for the density grid
    //! @param RESOLUTION the resolution to simulate for
    SimulateGaussianSphere::SimulateGaussianSphere( const linal::Vector3D &GRID_SPACING, const double RESOLUTION) :
      m_GridSpacing( GRID_SPACING),
      m_Resolution( RESOLUTION),
      m_Margin( 2)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SimulateGaussianSphere
    SimulateGaussianSphere *SimulateGaussianSphere::Clone() const
    {
      return new SimulateGaussianSphere( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SimulateGaussianSphere::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution
    //! @param RESOLUTION the resolution for the density map to be generated
    void SimulateGaussianSphere::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set the resolution
    double SimulateGaussianSphere::GetResolution() const
    {
      return m_Resolution;
    }

    //! @brief set the grid spacing
    //! @param GRID_SPACING the width of a grid element in x, y and z
    void SimulateGaussianSphere::SetGridSpacing( const linal::Vector3D &GRID_SPACING)
    {
      m_GridSpacing = GRID_SPACING;
    }

    //! @brief set the margin
    //! @param MARGIN number of additional cells next to last atom occupied cells
    void SimulateGaussianSphere::SetMargin( const size_t MARGIN)
    {
      m_Margin = MARGIN;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate simulated density from given list of atoms
    //! @param ATOMS siptrvector of atoms
    //! @return a simulated density map
    Map SimulateGaussianSphere::operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const
    {
      // measure ATOMS extent
      storage::VectorND< 2, linal::Vector3D> min_max_coord( DetermineGridCorners( ATOMS));

      linal::Vector3D &mincoord( min_max_coord.First());
      linal::Vector3D &maxcoord( min_max_coord.Second());

      // for undefined corners, return empty map
      if( !mincoord.IsDefined() || !maxcoord.IsDefined())
      {
        return Map();
      }

      // add margin
      mincoord -= m_Margin * m_GridSpacing;
      maxcoord += m_Margin * m_GridSpacing;

      // determine index
      const linal::VectorND< int, 3> index
      (
        int( std::floor( mincoord.X() / m_GridSpacing.X())),
        int( std::floor( mincoord.Y() / m_GridSpacing.Y())),
        int( std::floor( mincoord.Z() / m_GridSpacing.Z()))
      );

      // dimensions of grid
      const storage::VectorND< 3, size_t> dimensions
      (
        size_t( std::ceil( maxcoord.X() / m_GridSpacing.X())) - index( 0) ,
        size_t( std::ceil( maxcoord.Y() / m_GridSpacing.Y())) - index( 1),
        size_t( std::ceil( maxcoord.Z() / m_GridSpacing.Z())) - index( 2)
      );
      // allocate mask size and set values to 0
      math::Tensor< double> grid
      (
        dimensions.Third(), dimensions.Second(), dimensions.First(), double( 0.0)
      );

      // store the real space index for easier access
      const linal::Vector3D realspaceindex
        (
          index( 0)  * m_GridSpacing.X(),
          index( 1) * m_GridSpacing.Y(),
          index( 2)  * m_GridSpacing.Z()
        );

      // constants describing gaussian blob shape
      const double blob_k( math::Sqr( math::g_Pi / ( 2.4 + 0.8 * m_Resolution)));
      const double blob_c( math::Pow( blob_k / math::g_Pi, 1.5));

      // cutoff distance at 3 sigma
      const double cutoff_square( math::Sqr( 3.0 * ( 1.0 / math::Sqrt( 2.0)) * ( ( 2.4 + 0.8 * m_Resolution) / math::g_Pi)));
      const double cutoff( math::Sqrt( cutoff_square));

      // iterate over all atoms
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr
      )
      {
        const linal::Vector3D &current_coord( ( *itr)->GetCoordinates());

        // skip undefined coordinates
        if( !current_coord.IsDefined())
        {
          continue;
        }

        //  calculate the position within the grid
        const linal::Vector3D coord_in_grid( current_coord - realspaceindex);

        // mass of atom times c constant for sphere/blob
        const double blob_c_mass( blob_c * ( *itr)->GetType()->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass));

        // min and max index within grid
        const storage::VectorND< 3, size_t> min_index
        (
          size_t( std::max( int( 0), int( std::floor( ( coord_in_grid.X() - cutoff) / m_GridSpacing.X())))),
          size_t( std::max( int( 0), int( std::floor( ( coord_in_grid.Y() - cutoff) / m_GridSpacing.Y())))),
          size_t( std::max( int( 0), int( std::floor( ( coord_in_grid.Z() - cutoff) / m_GridSpacing.Z()))))
        );
        const storage::VectorND< 3, size_t> max_index
        (
          size_t( std::min( int( dimensions.First()), int( std::ceil( ( coord_in_grid.X() + cutoff) / m_GridSpacing.X())))),
          size_t( std::min( int( dimensions.Second()), int( std::ceil( ( coord_in_grid.Y() + cutoff) / m_GridSpacing.Y())))),
          size_t( std::min( int( dimensions.Third()), int( std::ceil( ( coord_in_grid.Z() + cutoff) / m_GridSpacing.Z()))))
        );

        // iterate over voxel in grid to calculate the intensity/density
        for( size_t i( min_index.First()); i < max_index.First(); ++i)
        {
          linal::Vector3D pos_voxel;
          pos_voxel.X() = i * m_GridSpacing.X();
          for( size_t j( min_index.Second()); j < max_index.Second(); ++j)
          {
            pos_voxel.Y() = j * m_GridSpacing.Y();
            for( size_t k( min_index.Third()); k < max_index.Third(); ++k)
            {
              pos_voxel.Z() = k * m_GridSpacing.Z();

              const double square_distance( ( coord_in_grid - pos_voxel).SquareNorm());

              // check distance cutoff
              if( square_distance <= cutoff_square)
              {
                // calculate the distance between point in mask grid and specified coordinate
                grid( k, j, i) += blob_c_mass * std::exp( -blob_k * square_distance);
              }
            } // x
          } // y
        } // z
      } // iterate over atoms

      // intervals
      linal::VectorND< int, 3> intervals;
      for( size_t i( 0); i < 3; ++i)
      {
        if( index( i) <= -int( index( i)))
        {
          intervals( i) = math::Absolute( index( i));
        }
        else if( index( i) >= 0)
        {
          intervals( i) = index( i) + index( i) - 1;
        }
        else
        {
          intervals( i) = index( i) - 1;
        }
      }

      // length
      const linal::Vector3D length
      (
        m_GridSpacing.X() * double( intervals( 0)),
        m_GridSpacing.Y() * double( intervals( 1)),
        m_GridSpacing.Z() * double( intervals( 2))
      );

      // end
      return Map
             (
               grid,
               index,
               intervals,
               length,
               m_GridSpacing,
               Map::GetDefaultAngle(),
               Map::GetDefaultAxis(),
               linal::Vector3D( 0.0)
             );
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SimulateGaussianSphere::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_GridSpacing    , ISTREAM);
      io::Serialize::Read( m_Resolution     , ISTREAM);
      io::Serialize::Read( m_Margin         , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SimulateGaussianSphere::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_GridSpacing    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Resolution     , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Margin         , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace density
} // namespace bcl
