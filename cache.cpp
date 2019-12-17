#include "cache.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;

// GLOBAL VARIABLES

cache_t *L1_Cache;
cache_t *L2_Cache;
cache_t *Victim_Cache;

vector<uint64_t> block_addresses;

int total_bits = 64;
struct cache_config_t *global_conf;

/** @brief Function to initialize your cache structures and any globals that you might need
 *
 *  @param conf pointer to the cache configuration structure
 *
 */
void cache_init(struct cache_config_t *conf)
{
    global_conf = conf;
    L1_Cache = new cache_t;
    L2_Cache = new cache_t;
    Victim_Cache = new cache_t;
    // GENERAL INFO NEEDED
    uint64_t blocksize_bits = conf->b;
    uint64_t bytes_per_block = 1<<blocksize_bits;

    // L1 SETUP

    L1_Cache->cachesize = 1<<(conf->c);
    uint64_t l1_amount_of_blocks = L1_Cache->cachesize/bytes_per_block;
    // uint64_t l1_blocks_per_set = 1<<(conf->s);
    L1_Cache->ways = 1<<(conf->s);

    int l1_amount_of_lines = ((L1_Cache->cachesize)/bytes_per_block)/(L1_Cache->ways);

    L1_Cache->offset_bits = conf->b;
    L1_Cache->index_bits = log2(l1_amount_of_lines);
    L1_Cache->tag_bits = total_bits - L1_Cache->offset_bits - L1_Cache->index_bits;

    // L1 INITIALIZATION

    L1_Cache->DataXLength = l1_amount_of_lines;
    L1_Cache->data = create2DCacheArray(l1_amount_of_lines, L1_Cache->ways);

    // L2 SETUP

    L2_Cache->cachesize = 1<<(conf->C);
    uint64_t l2_amount_of_blocks = L2_Cache->cachesize/bytes_per_block;
    L2_Cache->ways = 1<<(conf->S);

    int l2_amount_of_lines = (L2_Cache->cachesize/bytes_per_block)/L2_Cache->ways;

    L2_Cache->offset_bits = conf->b;
    L2_Cache->index_bits = log2(l2_amount_of_lines);
    L2_Cache->tag_bits = total_bits - L2_Cache->offset_bits - L2_Cache->index_bits;

    // L2 INITIALIZATION

    L2_Cache->DataXLength = l2_amount_of_lines;
    L2_Cache->data = create2DCacheArray(l2_amount_of_lines, L2_Cache->ways);

    // VICTIM SETUP

    Victim_Cache->cachesize = conf->v * bytes_per_block;
    Victim_Cache->ways = conf->v;
    Victim_Cache->offset_bits = conf->b;
    Victim_Cache->index_bits = 0; // log(1) = 0
    Victim_Cache->tag_bits = total_bits - Victim_Cache->offset_bits;

    // VICTIM INITIALIZATION

    Victim_Cache->DataXLength = 1;
    Victim_Cache->data = create2DCacheArray(1, Victim_Cache->ways);
    Victim_Cache->using_lru = 0;

}

/** @brief Function to initialize your cache structures and any globals that you might need
 *
 *  @param addr The address being accessed
 *  @param rw Tell if the access is a read or a write
 *  @param stats Pointer to the cache statistics structure
 *
 */
void cache_access(uint64_t addr, char rw, struct cache_stats_t *stats)
{
    int is_write = rw == WRITE;

    bool is_compulsory = true;

    /*

    // TESTS

    // IS COMPULSORY
    if (std::find(block_addresses.begin(), block_addresses.end(), getBlockAddress(addr, L1_Cache)) != block_addresses.end()) {
        is_compulsory = false;
    }

    // IS SIMILAR
    if (is_compulsory) {
        for (vector<uint64_t>::iterator it = block_addresses.begin() ; it != block_addresses.end(); ++it) {
            cout << *it;
            if ((*it - addr < 1<<20) || (addr - *it < 1<<20)) {
                cout << "couldhavebeenoptimized\n";
            }
        }
    }

    */

    // L1 LOOKUP

    uint64_t l1_index = getIndex(addr, L1_Cache);
    uint64_t l1_tag = getTag(addr, L1_Cache);

    // if (!(stats->num_accesses < 1000) && !is_compulsory) cout << "non comp miss\n";

    stats->num_accesses++;

    if (is_write) stats->num_accesses_writes++; else stats->num_accesses_reads++;

    int l1_possible_way = is_in_block(l1_index, l1_tag, addr, L1_Cache);
    if (l1_possible_way != -1) { // L1 HIT
        if (stats->num_accesses < 1000) cout << "hit\n";
        int way = l1_possible_way;

        if (L1_Cache->data[l1_index][way]->was_prefetched) {
            stats->num_useful_prefetches++;
            L1_Cache->data[l1_index][way]->was_prefetched = 0;
        } 

        // Write Block
        cache_write(addr, l1_index, l1_tag, way, is_write, L1_Cache, stats, 0, 0);
    } else { // L1 Miss!
        stats->num_misses_l1++;
        if (stats->num_accesses < 1000) {
            if (is_compulsory) {
                cout << "compulsory miss\n";
            } else {
                cout << "miss\n";
            }
        }
        if (is_write) stats->num_misses_writes_l1++; else stats->num_misses_reads_l1++;

        // VC LOOKUP

        uint64_t vc_tag = getTag(addr, Victim_Cache);
        int vc_possible_way = is_in_block(0, vc_tag, addr, Victim_Cache);

        if (vc_possible_way != -1) { // VC HIT
            stats->num_hits_vc++;
            int way = vc_possible_way;

            if (Victim_Cache->data[0][way]->was_prefetched) {
                stats->num_useful_prefetches++;
                Victim_Cache->data[0][way]->was_prefetched = 0;
            }

            //Victim_Cache->data[0][way]->valid = 0;
            // lru move to back
            LRUMoveToBack(0, way, Victim_Cache);

            int vc_dirty = Victim_Cache->data[0][way]->dirty;
            int should_write = is_write || vc_dirty;

            int l1_way = best_way_to_replace(l1_index, L1_Cache);

            // Writes Block to L1 Cache
            cache_write(addr, l1_index, l1_tag, l1_way, should_write, L1_Cache, stats, 0, 0);

            // before or after?
        } else { // VC MISS
            stats->num_misses_vc++;
            if (is_write) stats->num_misses_writes_vc++; else stats->num_misses_reads_vc++;

            // L2 LOOKUP

            uint64_t l2_index = getIndex(addr, L2_Cache);
            uint64_t l2_tag = getTag(addr, L2_Cache);

            int l2_possible_way = is_in_block(l2_index, l2_tag, addr, L2_Cache);
            if (l2_possible_way != -1) { // L2 HIT
                int way = l2_possible_way;

                if (L2_Cache->data[l2_index][way]->was_prefetched) {
                    stats->num_useful_prefetches++;
                    L2_Cache->data[l2_index][way]->was_prefetched = 0;
                }

                // cache_write(addr, l2_index, l2_tag, way, 0, L2_Cache, stats, 0);
                // instead update LRU
                LRUMoveToFront(l2_index, way, L2_Cache);

                int l2_dirty = L2_Cache->data[l2_index][way]->dirty;
                int should_write = is_write || l2_dirty;

                // not sure about this line
                L2_Cache->data[l2_index][way]->dirty = 0;

                int l1_way = best_way_to_replace(l1_index, L1_Cache);
                cache_write(addr, l1_index, l1_tag, l1_way, should_write, L1_Cache, stats, 0, 0); // IS_WRITE??? or 0?
                
            } else { // L2 MISS
                stats->num_misses_l2++;
                if (is_write) stats->num_misses_writes_l2++; else stats->num_misses_reads_l2++;

                // must get from memory!

                int l2_way = best_way_to_replace(l2_index, L2_Cache);
                cache_write(addr, l2_index, l2_tag, l2_way, 0, L2_Cache, stats, 0, 0);

                stats->num_bytes_transferred += 1<<(L2_Cache->offset_bits);

                int l1_way = best_way_to_replace(l1_index, L1_Cache);
                cache_write(addr, l1_index, l1_tag, l1_way, is_write, L1_Cache, stats, 0, 0);

                // PREFETCH

                prefetch2(addr, stats); 
            }
        }
        
    }
    if (is_compulsory) block_addresses.push_back(getBlockAddress(addr, L1_Cache));
}

/** @brief Function to free any allocated memory and finalize statistics
 *
 *  @param stats pointer to the cache statistics structure
 *
 */
void cache_cleanup(struct cache_stats_t *stats)
{
    free_memory();
    stats->miss_rate_l1 = (double)stats->num_misses_l1/(double)stats->num_accesses;
    stats->miss_rate_l2 = (double)stats->num_misses_l2/(double)stats->num_misses_vc;
    stats->miss_rate_vc = (double)stats->num_misses_vc/(double)stats->num_misses_l1;
    double miss_penalty_l1_vc = (stats->hit_time_l2) + (stats->miss_rate_l2 * stats->hit_time_mem);
    double avg_access_time_l1 = (stats->hit_time_l1) + (stats->miss_rate_l1 * stats->miss_rate_vc)*(miss_penalty_l1_vc);
    stats->avg_access_time = avg_access_time_l1;
}

void cache_write(uint64_t addr, uint64_t index, uint64_t tag, int way, int is_write, cache_t *cache, struct cache_stats_t *stats, int ignore_lru, int is_prefetch) {
    cache_entity_t * cache_entity = cache->data[index][way];
    // clean up/replacement
    if (cache == L1_Cache) {
        LRUMoveToFront(index, way, L1_Cache);
        if (cache_entity->valid && cache_entity->tag != tag) { // push to VC/L2
            if (global_conf->v != 0) { // if VC exists, push to vc
                uint64_t vc_tag = getTag(cache_entity->address, Victim_Cache);
                if (is_in_block(0, vc_tag, cache_entity->address, Victim_Cache) == -1) {
                    int vc_way = best_way_to_replace(0, Victim_Cache);                
                    cache_write(cache_entity->address, 0, vc_tag, vc_way, cache_entity->dirty, Victim_Cache, stats, 0, 0);
                } else {
                    printf("In block already??");
                }
            } else { // if vc does not exist, push to L2 but only if dirty
                if (cache_entity->dirty) {
                    // write to L2
                    uint64_t l2_index = getIndex(cache_entity->address, L2_Cache);
                    uint64_t l2_tag = getTag(cache_entity->address, L2_Cache);
                    int l2_possible_way = is_in_block(l2_index, l2_tag, cache_entity->address, L2_Cache);
                    if (l2_possible_way != -1) { // if already exists just update dirty bit
                        L2_Cache->data[l2_index][l2_possible_way]->dirty = 1;
                        //LRUMoveToFront(l2_index, l2_possible_way, L2_Cache);
                    } else {
                        int way_to_replace = best_way_to_replace(l2_index, L2_Cache);
                        cache_write(cache_entity->address, l2_index, l2_tag, way_to_replace, 1, L2_Cache, stats, 1, 0);
                        //LRUMoveToFront(l2_index, way_to_replace, L2_Cache);
                    }
                }
            }
        }
        if (cache_entity->tag != tag && cache_entity->dirty) cache_entity->dirty = 0; // set to 0 so it can be set differently
    } else if (cache == Victim_Cache) {
        if (cache_entity->valid && cache_entity->tag != tag) {
            LRUMoveToFront(index, way, Victim_Cache);
            if (cache_entity->dirty) {
                // write to L2
                uint64_t l2_index = getIndex(cache_entity->address, L2_Cache);
                uint64_t l2_tag = getTag(cache_entity->address, L2_Cache);
                int l2_possible_way = is_in_block(l2_index, l2_tag, cache_entity->address, L2_Cache);
                if (l2_possible_way != -1) { // if already exists just update dirty bit
                    L2_Cache->data[l2_index][l2_possible_way]->dirty = 1;
                    //LRUMoveToFront(l2_index, l2_possible_way, L2_Cache);
                } else {
                    int way_to_replace = best_way_to_replace(l2_index, L2_Cache);
                    cache_write(cache_entity->address, l2_index, l2_tag, way_to_replace, 1, L2_Cache, stats, 1, 0);
                    //LRUMoveToFront(l2_index, way_to_replace, L2_Cache);
                }
            }
        } else if (cache_entity->valid && cache_entity->tag == tag) {
            printf("In block already?????");
        }
        if (cache_entity->tag != tag && cache_entity->dirty) cache_entity->dirty = 0; // set to 0 so it can be set differently
    } else if (cache == L2_Cache) {
        if (ignore_lru == 0) LRUMoveToFront(index, way, L2_Cache);

        //if (is_prefetch) LRUMoveToBack(index, way, L2_Cache);
        
        if (cache_entity->valid && cache_entity->tag != tag) { // write back
            if (cache_entity->dirty) {
                stats->num_write_backs++;
                stats->num_bytes_transferred += 1<<(cache->offset_bits);
            } 
        }
        if (cache_entity->tag != tag && cache_entity->dirty) cache_entity->dirty = 0; // set to 0 so it can be set differently
    } 

    // Do ready/write op

    if ((cache_entity->tag == tag && cache_entity->dirty) || is_write) { // sets dirty bit if necessary
        cache_entity->dirty = 1;
    } else {
        cache_entity->dirty = 0;
    }

    cache_entity->valid = 1; // sets valid bit to true
    cache_entity->tag = tag;
    cache_entity->address = addr;
    cache_entity->was_prefetched = is_prefetch;
}

void prefetch(uint64_t addr, struct cache_stats_t *stats) {
    uint64_t block_address = getBlockAddress(addr, L2_Cache);
    for (uint64_t i = 1; i <= global_conf->k; i++) {
        uint64_t block_to_prefetch = (block_address+i)<<global_conf->b;
        uint64_t index_to_prefetch = getIndex(block_to_prefetch, L2_Cache);
        uint64_t tag_to_prefetch = getTag(block_to_prefetch, L2_Cache);
        int best_way = best_way_to_replace(index_to_prefetch, L2_Cache);
        int is_in_the_block = is_in_block(index_to_prefetch, tag_to_prefetch, addr, L2_Cache);
        if (is_in_the_block == -1) {
            cache_write(block_to_prefetch, index_to_prefetch, tag_to_prefetch, best_way, 0, L2_Cache, stats, 1, 1);
            L2_Cache->data[index_to_prefetch][best_way]->was_prefetched = 1;
            stats->num_prefetches++;
            stats->num_bytes_transferred += 1<<(L2_Cache->offset_bits);
        }
    }
}

void prefetch2(uint64_t addr, struct cache_stats_t *stats) {
    uint64_t difference = (1<<(L2_Cache->offset_bits));
    for (uint64_t block_to_prefetch = addr + difference; block_to_prefetch <= addr + (global_conf->k)*difference; block_to_prefetch += difference) {
        uint64_t index_to_prefetch = getIndex(block_to_prefetch, L2_Cache);
        uint64_t tag_to_prefetch = getTag(block_to_prefetch, L2_Cache);
        int best_way = best_way_to_replace(index_to_prefetch, L2_Cache);
        int is_in_the_block = is_in_block(index_to_prefetch, tag_to_prefetch, addr, L2_Cache);
        if (is_in_the_block == -1) {
            cache_write(block_to_prefetch, index_to_prefetch, tag_to_prefetch, best_way, 0, L2_Cache, stats, 1, 1);
            L2_Cache->data[index_to_prefetch][best_way]->was_prefetched = 1;
            stats->num_prefetches++;
            stats->num_bytes_transferred += 1<<(L2_Cache->offset_bits);
        }
    }
}

int best_way_to_replace(uint64_t index, cache_t *cache) {
    int best_way;
    for (int i = 0; i < cache->ways; i++) {
        if (cache->data[index][i]->valid == 0) {
            best_way = i;
            return best_way;
        } else {
            if (cache->data[index][i]->priority == cache->ways - 1) { // has worst priority
                best_way = i;
                return best_way;
            }
        }
    }
}

int best_way_to_replace2(uint64_t index, cache_t *cache) {
    int best_way;
    for (int i = 0; i < cache->ways; i++) {
        if (cache->data[index][i]->valid == 0) {
            best_way = i;
            return best_way;
        }
    }

    for (int i = 0; i < cache->ways; i++) {
        if (cache->data[index][i]->priority == cache->ways - 1) { // has worst priority
            best_way = i;
            return best_way;
        }
    }
}

// Get bits tag
uint64_t getTag(uint64_t address, cache_t *cache) {
  // DEBUG: printf("\nAddress: %x", address);
  return (getBits(address, cache->offset_bits + cache->index_bits, cache->tag_bits));
}

// Get bits index
uint64_t getIndex(uint64_t address, cache_t *cache) {
  uint64_t index = getBits(address, cache->offset_bits, cache->index_bits);
  return index;
}

// Get block address
uint64_t getBlockAddress(uint64_t address, cache_t *cache) {
    uint64_t block_address = getBits(address, cache->offset_bits, cache->index_bits + cache->tag_bits);
    return block_address;
}

// returns the way index if the tag is located in the specified cache. else returns 0
int is_in_block(uint64_t index, uint64_t tag, uint64_t addr, cache_t *cache) {
  for (int way = 0; way < cache->ways; way++) {
    if ((cache->data)[index][way]->tag == tag && (cache->data)[index][way]->valid) {
      return way;
    }
  }
  return -1;
}

// moves the way to the front of the LRU (if it is read or first written) 1 (1) 2 (2) 3 (3) 4 
int LRUMoveToFront(unsigned index, int way, cache_t *cache) {
  // move around the priority
  int ways = cache->ways;
  int old_priority = cache->data[index][way]->priority;
  for (int i = 0; i < ways; i++) {
      if (i == way) {
          (cache->data)[index][i]->priority = 0;
      } else if (cache->data[index][i]->priority < old_priority) { // increments priority only if it had better priority than the way (lower number)
          (cache->data)[index][i]->priority++;
      }
  }
}

// moves the way to the back of the LRU. Maybe doesn't work because index 0 is invalid, 1,2,3 are valid, 4, 5, 6 not. index 0 gets replaced, but is placed super super behind
int LRUMoveToBack(unsigned index, int way, cache_t *cache) {
  // move around the priority
  int ways = cache->ways;
  int old_priority = cache->data[index][way]->priority;
  for (int i = 0; i < ways; i++) {
      if (i == way) {
          (cache->data)[index][i]->priority = cache->ways - 1;
      } else if (cache->data[index][i]->priority > old_priority) { // decrements priority only if it had better priority than the way (lower number)
          (cache->data)[index][i]->priority--;
      }
  }
}

int FIFOSendBack(uint64_t index, int way, cache_t *cache) { // 1 (2) 2 (3) 3 (4) 4 (1)
    int ways = cache->ways;
    for (int i = 0; i < ways; i++) {
        if (i == way) {
            cache->data[index][i]->priority = way;
        } else if (i > way) {
            cache->data[index][i]->priority--;
        }
    }
}

int moveToFrontOfArray(unsigned number, unsigned length, int *Array) {
  // get index of number
  int num_index;
  for (int i = 0; i < length; i++) {
    if (Array[i] == number) num_index = i;
  }

  if (num_index == 0) return 1;

  for (int num_to_move = num_index-1; num_to_move >= 0; num_to_move--) {
    Array[num_to_move + 1] = Array[num_to_move];
  }
  Array[0] = number;
}

// Helper functions

unsigned createMask(unsigned a, unsigned b)
{
   unsigned r = 0;
   for (unsigned i=a; i<=b; i++)
       r |= 1 << i;

   return r;
}

uint64_t getBits( uint64_t allbits, unsigned lsb, unsigned count )
{
    uint64_t mask = ~(0xffffffffffffffff << count) << lsb ;
    return (allbits & mask) >> lsb ;
}

uint64_t ***create3DIndexArray(int x, int y, int z) {
  uint64_t*** A = (uint64_t***)malloc(x * sizeof(uint64_t**));

	if (A == NULL) {
		fprintf(stderr, "\nOut of memory");
		exit(0);
	}

	for (int i = 0; i < x; i++)
	{
		A[i] = (uint64_t**)malloc(y * sizeof(uint64_t*));
		if (A[i] == NULL) {
			fprintf(stderr, "\nOut of memory");
			exit(0);
		}

		for (int j = 0; j < y; j++)
		{
			A[i][j] = (uint64_t*)malloc(z * sizeof(uint64_t));
	   		if (A[i][j] == NULL) {
				fprintf(stderr, "\nOut of memory");
				exit(0);
			}
		}
	}

  // assign NULL values to allocated memory
	for (int i = 0; i < x; i++)
    for (int j = 0; j < y; j++)
      for (int k = 0; k < z; k++)
        A[i][j][k] = 0;

  // assign values to allocated memory
	for (int i = 0; i < x; i++)
    for (int k = 0; k < z; k++)
      A[i][1][k] = 0; // sets the valid bit to 0 on initialization

  return A;
}

int **create2DLRUArray(int x, int y) {
  int**A = (int**)malloc(x * sizeof(int*));

	if (A == NULL) {
		fprintf(stderr, "\nOut of memory");
		exit(0);
	}

	for (int i = 0; i < x; i++)
	{
		A[i] = (int*)malloc(y * sizeof(int*));
		if (A[i] == NULL) {
			fprintf(stderr, "\nOut of memory");
			exit(0);
		}
	}

  // assign values to allocated memory
	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
      A[i][j] = j; // sets the valid bit to 0 on initialization

  return A;
}

cache_entity_t ***create2DCacheArray(int x, int y) {
  cache_entity_t***A = (cache_entity_t***)malloc(x * sizeof(cache_entity_t**));

	if (A == NULL) {
		fprintf(stderr, "\nOut of memory");
		exit(0);
	}

	for (int i = 0; i < x; i++)
	{
		A[i] = (cache_entity_t**)malloc(y * sizeof(cache_entity_t**));
		if (A[i] == NULL) {
			fprintf(stderr, "\nOut of memory");
			exit(0);
		}
	}

  // assign values to allocated memory
	for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            cache_entity_t * cache_entity = new cache_entity_t;
            cache_entity->valid = 0;
            cache_entity->priority = j;
            cache_entity->was_prefetched = 0;
            cache_entity->tag = NULL;
            cache_entity->address = NULL;
            A[i][j] = cache_entity;
        }
    }

  return A;
}

void free_memory() {
    for (int i = 0; i < L1_Cache->DataXLength; i++) {
        free(L1_Cache->data[i]);
    }

    for (int i = 0; i < Victim_Cache->DataXLength; i++) {
        free(Victim_Cache->data[i]);
    }

    for (int i = 0; i < L2_Cache->DataXLength; i++) {
        free(L2_Cache->data[i]);
    }

    free(L1_Cache->data);
    free(Victim_Cache->data);
    free(L2_Cache->data);
}