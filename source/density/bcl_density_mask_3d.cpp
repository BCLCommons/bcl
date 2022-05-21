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
#include "density/bcl_density_mask_3d.h"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_map.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Mask3d::s_Instance
    (
      GetObjectInstances().AddInstance( new Mask3d())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Mask3d::Mask3d() :
      m_Mask(),
      m_Index(),
      m_Position(),
      m_GridSpacing()
    {
    }

    //! @brief construct from list of coordinates
    //! @param COORDS a list of coordinates that defines the mask
    //! @param MASKING_DISTANCE distance for sigmoid in angstrom function to have marginal impact
    //! @param GRID_SPACING length of cube in grid
    //! @param CELL_POSITION position of the referenced grid in space
    Mask3d::Mask3d
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDS,
      const double MASKING_DISTANCE,
      const linal::Vector3D &GRID_SPACING,
      const linal::Vector3D &CELL_POSITION
    ) :
      m_Mask( 0, 0, 0),
      m_Index( 0, 0, 0),
      m_Position( CELL_POSITION),
      m_GridSpacing( GRID_SPACING)
    {
      // determine the corners as defined by given coordinates
      storage::VectorND< 2, linal::Vector3D> corners( DetermineGridCorners( COORDS));

      // undefined corners
      if( !corners.First().IsDefined())
      {
        BCL_MessageDbg( "could not determine corners for mask");
        return;
      }

      // add margin
      corners = AddMargin( corners, MASKING_DISTANCE);

      // find index of mask relative to map
      m_Index.First()  = int( std::floor( ( corners.First().X() - m_Position.X()) / GRID_SPACING.X()));
      m_Index.Second() = int( std::floor( ( corners.First().Y() - m_Position.Y()) / GRID_SPACING.Y()));
      m_Index.Third()  = int( std::floor( ( corners.First().Z() - m_Position.Z()) / GRID_SPACING.Z()));

      // allocate mask size and set values to 1
      m_Mask = math::Tensor< double>
      (
        size_t( std::ceil( ( corners.Second().X() - m_Position.X()) / GRID_SPACING.X())) - m_Index.First() ,
        size_t( std::ceil( ( corners.Second().Y() - m_Position.Y()) / GRID_SPACING.Y())) - m_Index.Second(),
        size_t( std::ceil( ( corners.Second().Z() - m_Position.Z()) / GRID_SPACING.Z())) - m_Index.Third() ,
        double( 1.0)
      );

      // store the real space index for easier access
      const linal::Vector3D realspaceindex
        (
          m_Position.X() + m_Index.First() * m_GridSpacing.X(),
          m_Position.Y() + m_Index.Second() * m_GridSpacing.Y(),
          m_Position.Z() + m_Index.Third() * m_GridSpacing.Z()
        );

      linal::Vector3D pos_voxel;

      // iterate over every point in grid to calculate the masking value
      for( size_t i( 0); i < m_Mask.GetNumberCols(); ++i)
      {
        pos_voxel.X() = i * m_GridSpacing.X() + realspaceindex.X();
        for( size_t j( 0); j < m_Mask.GetNumberRows(); ++j)
        {
          pos_voxel.Y() = j * m_GridSpacing.Y() + realspaceindex.Y();
          for( size_t k( 0); k < m_Mask.NumberLayers(); ++k)
          {
            pos_voxel.Z() = k * m_GridSpacing.Z() + realspaceindex.Z();

            double product( 1.0);
            // iterate over all coordinates
            for
            (
              util::SiPtrVector< const linal::Vector3D>::const_iterator itr( COORDS.Begin()), itr_end( COORDS.End());
              itr != itr_end;
              ++itr
            )
            {
              const linal::Vector3D &current_coord( **itr);
              if( !current_coord.IsDefined())
              {
                continue;
              }

              // calculate the distance between point in mask grid and specified coordinate
              product *= 1 - Sigmoid( MASKING_DISTANCE - ( current_coord - pos_voxel).Norm());
            }

            // actual mask weight is 1 minus the calculated product
            m_Mask( k, j, i) -= product;
          }
        }
      }
    };

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the cross correlation between two density maps weighted by mask
    //! @param DENSITY_MAP_EXP experimental density map
    //! @param DENSITY_MAP_SIM simulated density map
    //! @return CCC
    double Mask3d::CrossCorrelationCoefficient
    (
      const Map &DENSITY_MAP_EXP,
      const Map &DENSITY_MAP_SIM
    ) const
    {
      // assert the all maps' cell widths are identical
      BCL_Assert
      (
           DENSITY_MAP_EXP.GetCellWidth() == DENSITY_MAP_SIM.GetCellWidth()
        && DENSITY_MAP_EXP.GetCellWidth() == m_GridSpacing,
        "The density maps and the mask need to have the same Grid - spacing:/nmask: "
        + util::Format()(m_GridSpacing)
        + "\nsimulated density map: "
        + util::Format()( DENSITY_MAP_SIM.GetCellWidth())
        + "\nexperimental density map: "
        + util::Format()( DENSITY_MAP_EXP.GetCellWidth())
      );

      // determine overlap of each density map with this mask
      const storage::VectorND< 3, int> dm_rel_index_sim( RelativeIndex( DENSITY_MAP_SIM));
      const storage::VectorND< 3, int> dm_rel_index_exp( RelativeIndex( DENSITY_MAP_EXP));

      // determine the common overlap of the two densities and this mask
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > commonarea
      (
        CommonOverlap( OverlappingIndices( DENSITY_MAP_EXP), OverlappingIndices( DENSITY_MAP_SIM))
      );

      // mean and sd of experimental and simulated map
      const math::RunningAverageSD< double> meansd_exp( CalculateMeanSDCommonRegion( DENSITY_MAP_EXP, commonarea));
      const math::RunningAverageSD< double> meansd_sim( CalculateMeanSDCommonRegion( DENSITY_MAP_SIM, commonarea));
      const double mean_exp( meansd_exp.GetAverage());
      const double mean_sim( meansd_sim.GetAverage());
      double coefficient( 0);

      // iterate over mask dimensions
      for( size_t i( commonarea.First().First()); i < commonarea.Second().First(); ++i)
      {
        const size_t exp_index_x( i - dm_rel_index_exp.First());
        const size_t sim_index_x( i - dm_rel_index_sim.First());
        for( size_t j( commonarea.First().Second()); j < commonarea.Second().Second(); ++j)
        {
          const size_t exp_index_y( j - dm_rel_index_exp.Second());
          const size_t sim_index_y( j - dm_rel_index_sim.Second());
          for( size_t k( commonarea.First().Third()); k < commonarea.Second().Third(); ++k)
          {
            const size_t exp_index_z( k - dm_rel_index_exp.Third());
            const size_t sim_index_z( k - dm_rel_index_sim.Third());
            // sum up products
            coefficient += m_Mask( k, j, i)
              *
              (
                DENSITY_MAP_EXP
                (
                  exp_index_x,
                  exp_index_y,
                  exp_index_z
                ) - mean_exp
              )
              *
              (
                DENSITY_MAP_SIM
                (
                  sim_index_x,
                  sim_index_y,
                  sim_index_z
                ) - mean_sim
              );
          }
        }
      }

      // number of voxels considered
      const size_t voxel_count
                (
                    ( commonarea.Second().First() - commonarea.First().First())
                  * ( commonarea.Second().Second() - commonarea.First().Second())
                  * ( commonarea.Second().Third() - commonarea.First().Third())
                );

      // calculate final correlation and return
      return coefficient * ( double( 1.0) / ( meansd_exp.GetStandardDeviation() * meansd_sim.GetStandardDeviation())) /
             voxel_count;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Mask3d::Read( std::istream &ISTREAM)
    {
      // read Index, Position and GridSpacing of the mask
      io::Serialize::Read( m_Index,        ISTREAM);
      io::Serialize::Read( m_Position,     ISTREAM);
      io::Serialize::Read( m_GridSpacing,  ISTREAM);

      // read the actual mask
      io::Serialize::Read( m_Mask,         ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Mask3d::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write Index, Position and GridSpacing of the mask
      io::Serialize::Write( m_Index,        OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Position,     OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GridSpacing,  OSTREAM, INDENT) << '\n';

      // write the actual mask
      io::Serialize::Write( m_Mask,         OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate mean and standard deviation of over the mask given a density map
    //! is not using the values in the mask, just the box that is defined by the mask
    //! @param DENSITY_MAP density map which's mean and SD shall be calculated for the area covered by the mask
    //! @return dataset statistics mean sd
    math::RunningAverageSD< double> Mask3d::CalculateMeanSD( const Map &DENSITY_MAP) const
    {
      // create empty data set statistic
      math::RunningAverageSD< double> mean_sd;

      // indices range
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > beg_end_indices( OverlappingIndices( DENSITY_MAP));

      // relative index of density map
      const storage::VectorND< 3, int> dm_rel_index( RelativeIndex( DENSITY_MAP));

      // iterate over given map over mask boundaries
      for( size_t i( beg_end_indices.First().First()); i < ( beg_end_indices.Second().First()); ++i)
      {
        const size_t index_x( i - dm_rel_index.First());
        for( size_t j( beg_end_indices.First().Second()); j < beg_end_indices.Second().Second(); ++j)
        {
          const size_t index_y( j - dm_rel_index.Second());
          for( size_t k( beg_end_indices.First().Third()); k < beg_end_indices.Second().Third(); ++k)
          {
            const size_t index_z( k - dm_rel_index.Third());
            // fill mean_sd with values
            mean_sd += DENSITY_MAP( index_x, index_y, index_z);
          }
        }
      }

      // end
      return mean_sd;
    }

    //! @brief determine corners of grid the coordinates fit in
    //! @param COORDS a list of coordinates that defines the grid
    //! @return min coord xyz and max coord xyz
    storage::VectorND< 2, linal::Vector3D> Mask3d::DetermineGridCorners
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDS
    )
    {
      // Initialize corners
      const linal::Vector3D nullcoord( std::numeric_limits< double>::max());
      storage::VectorND< 2, linal::Vector3D> corners( nullcoord, -nullcoord);

      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator itr( COORDS.Begin()), itr_end( COORDS.End());
        itr != itr_end;
        ++itr
      )
      {
        const linal::Vector3D &current_coord( **itr);
        if( !current_coord.IsDefined())
        {
          // Make sure, that the current coordinates are defined
          BCL_MessageDbg
          (
            "cannot use undefined coordinates" + util::Format()( current_coord)
          );
          continue;
        }

        // iterate over all atoms to find the maximum and minimum values
        corners.First().X() = std::min( current_coord.X(), corners.First().X());
        corners.First().Y() = std::min( current_coord.Y(), corners.First().Y());
        corners.First().Z() = std::min( current_coord.Z(), corners.First().Z());

        corners.Second().X() = std::max( current_coord.X(), corners.Second().X());
        corners.Second().Y() = std::max( current_coord.Y(), corners.Second().Y());
        corners.Second().Z() = std::max( current_coord.Z(), corners.Second().Z());
      }

      // no coord was defined
      if( corners.First().X() > corners.Second().X())
      {
        corners.First() = linal::Vector3D( util::GetUndefined< double>());
        corners.Second() = linal::Vector3D( util::GetUndefined< double>());
      }

      return corners;
    }

    //! @brief add margins to corners
    //! @param CORNERS min coord xyz and max coord xyz
    //! @param MARGIN margin to add
    //! @return corners with additional margin
    storage::VectorND< 2, linal::Vector3D> Mask3d::AddMargin
    (
      const storage::VectorND< 2, linal::Vector3D> &CORNERS,
      const double MARGIN
    )
    {
      storage::VectorND< 2, linal::Vector3D> result( CORNERS); //Variable where the result shall be saved in
      result.First() -= MARGIN;
      result.Second() += MARGIN;

      return result;
    }

    //! @brief determine the indices ranges that overlap with given density map, relative to this mask
    //! @param DENSITY_MAP overlapping density map
    //! @return pair of start indices and end indices
    storage::VectorND< 2, storage::VectorND< 3, size_t> > Mask3d::OverlappingIndices( const Map &DENSITY_MAP) const
    {
      // relative index of density map
      const storage::VectorND< 3, int> dm_rel_index( RelativeIndex( DENSITY_MAP));

      // dimenesion of density map
      const storage::VectorND< 3, size_t> dm_dim( DENSITY_MAP.GetDimensions());

      // find start index iteration
      const int x_beg( std::max( dm_rel_index.First() , int( 0)));
      const int y_beg( std::max( dm_rel_index.Second(), int( 0)));
      const int z_beg( std::max( dm_rel_index.Third() , int( 0)));

      // find maximum values for each coordinate
      const int x_end( std::min( dm_rel_index.First()  + int( dm_dim.First()) , int( m_Mask.GetNumberCols())));
      const int y_end( std::min( dm_rel_index.Second() + int( dm_dim.Second()), int( m_Mask.GetNumberRows())));
      const int z_end( std::min( dm_rel_index.Third()  + int( dm_dim.Third()) , int( m_Mask.NumberLayers())));

      // if mask does not overlap with given density map
      if( x_beg >= x_end || y_beg >= y_end || z_beg >= z_end)
      {
        return storage::VectorND< 2, storage::VectorND< 3, size_t> >
               (
                 storage::VectorND< 3, size_t>( 0, 0, 0), storage::VectorND< 3, size_t>( 0, 0, 0)
               );
      }
      else
      {
        return storage::VectorND< 2, storage::VectorND< 3, size_t> >
               (
                 storage::VectorND< 3, size_t>( size_t( x_beg), size_t( y_beg), size_t( z_beg)),
                 storage::VectorND< 3, size_t>( size_t( x_end), size_t( y_end), size_t( z_end))
               );
      }
    }

    //! @brief calculate relative index of given density map to the mask
    //! @param DENSITY_MAP in question
    //! @return relative index
    storage::VectorND< 3, int> Mask3d::RelativeIndex( const Map &DENSITY_MAP) const
    {
      // index of density map
      // return
      return storage::VectorND< 3, int>
             (
               DENSITY_MAP.GetIndex()( 0) - m_Index.First(),
               DENSITY_MAP.GetIndex()( 1) - m_Index.Second(),
               DENSITY_MAP.GetIndex()( 2) - m_Index.Third()
             );
    }

    math::RunningAverageSD< double> Mask3d::CalculateMeanSDCommonRegion
    (
      const Map &DENSITY_MAP,
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > &OVERLAP
    ) const
    {
      // create empty dataset
      math::RunningAverageSD< double> mean_sd;

      // relative index of density map
      const storage::VectorND< 3, int> dm_rel_index( RelativeIndex( DENSITY_MAP));

      for( size_t i( OVERLAP.First().First()); i < ( OVERLAP.Second().First()); ++i)
      {
        const size_t pos_x( i - dm_rel_index.First());
        for( size_t j( OVERLAP.First().Second()); j < OVERLAP.Second().Second(); ++j)
        {
          const size_t pos_y( j - dm_rel_index.Second());
          for( size_t k( OVERLAP.First().Third()); k < OVERLAP.Second().Third(); ++k)
          {
            const size_t pos_z( k - dm_rel_index.Third());
            // fill mean_sd with values
            mean_sd += DENSITY_MAP( pos_x, pos_y, pos_z);
          }
        }
      }
      return mean_sd;
    }

    //! @brief determines the overlap between the mask and two given ranges
    //! @param RANGE_ONE
    //! @param RANGE_TWO
    //! @return returns the common overlap
    storage::VectorND< 2, storage::VectorND< 3, size_t> > Mask3d::CommonOverlap
    (
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > &RANGE_ONE,
      const storage::VectorND< 2, storage::VectorND< 3, size_t> > &RANGE_TWO
    )
    {
      storage::VectorND< 2, storage::VectorND< 3, size_t> > overlap;
      overlap.First().First()  = std::max( RANGE_ONE.First().First() , RANGE_TWO.First().First()) ;
      overlap.First().Second() = std::max( RANGE_ONE.First().Second(), RANGE_TWO.First().Second());
      overlap.First().Third()  = std::max( RANGE_ONE.First().Third() , RANGE_TWO.First().Third()) ;

      overlap.Second().First()  = std::min( RANGE_ONE.Second().First() , RANGE_TWO.Second().First()) ;
      overlap.Second().Second() = std::min( RANGE_ONE.Second().Second(), RANGE_TWO.Second().Second());
      overlap.Second().Third()  = std::min( RANGE_ONE.Second().Third() , RANGE_TWO.Second().Third()) ;

      return overlap;
    }

  } // namespace density
} // namespace bcl
