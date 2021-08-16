/*
 * Copyright (C) 2015 University of Chicago.
 * See COPYRIGHT notice in top-level directory.
 *
 */

/*
* The test program generates some synthetic traffic patterns for the model-net network models.
* currently it only support the fat tree network model uniform random and nearest neighbor traffic patterns.
*/

#include "codes/model-net.h"
#include "codes/lp-io.h"
#include "codes/net/fattree.h"
#include "codes/codes.h"
#include "codes/codes_mapping.h"
#include "codes/configuration.h"
#include "codes/lp-type-lookup.h"

#define DEFAULT_PAYLOAD_SZ (512*256)

#define PARAMS_LOG 1

static int net_id = 0;
static int offset = 2;
static int traffic = 1;
static int payload_size = DEFAULT_PAYLOAD_SZ;
static double arrival_time = 1000.0;
static double load = 0.0;	//Percent utilization of terminal uplink
static double MEAN_INTERVAL = 0.0;
char * modelnet_stats_dir;
/* whether to pull instead of push */

static int num_servers_per_rep = 0;
static int num_routers_per_grp = 0;
static int num_nodes_per_grp = 0;

static int num_groups = 0;
static int num_nodes = 0;
static int num_msgs = 20;

typedef struct svr_msg svr_msg;
typedef struct svr_state svr_state;

/* global variables for codes mapping */
static char group_name[MAX_NAME_LENGTH];
static char lp_type_name[MAX_NAME_LENGTH];
static int group_index, lp_type_index, rep_id, offset;


/* convert GiB/s and bytes to ns */
 static tw_stime bytes_to_ns(uint64_t bytes, double GB_p_s)
 {
     tw_stime time;
 
     /* bytes to GB */
     time = ((double)bytes)/(1024.0*1024.0*1024.0);
     /* GiB to s */
     time = time / GB_p_s;
     /* s to ns */
     time = time * 1000.0 * 1000.0 * 1000.0;
 
     return(time);
 }

/* type of events */
enum svr_event
{
    KICKOFF,	   /* kickoff event */
    REMOTE,        /* remote event */
    LOCAL,      /* local event */
    ACK         /* event sent when remote endpoint has received message */
};
/* type of synthetic traffic */
enum TRAFFIC
{
	UNIFORM = 1, /* sends message to a randomly selected node */
    BISECTION = 2, /* sends messages to node established in bisection pairing*/
        RANK_ZERO = 3, /* Send all traffic to rank 0 */
        PFISTER = 4, /* A combination of UNIFORM and RANK_ZERO */
        GRID = 5, /* breaks into a cartesian grid and sends to 4 neighbors */
        UNIFORM_NN = 6,
        CSTENCIL = 7

};
enum TRAFFIC_MAP {
    REGULAR = 0,
    TILED = 1, 
    RANDOM = 2,
    RCM = 3,
    METIS = 4,
    IMPROVEMENT_PARTITION = 5,
    MAP_SCALE = 1000
};

struct svr_state
{
    int msg_sent_count;   /* requests sent */
    int msg_recvd_count;  /* requests recvd */
    int local_recvd_count; /* number of local messages received */
    tw_stime start_ts;    /* time that we started sending requests */
    tw_stime end_ts;      /* time that we ended sending requests */
    int recv_per_round;   /* In GRID, expected messages received per round */ 
};

    int hardcoded_mappings[6][4608] = {
        {}, {}, {}, 
        {
#include "/g/g14/taffet2/graph/rcm.csv"
        },
        {
#include "/g/g14/taffet2/graph/gpmetis-rb.csv"
        }, 
        {
#include "/g/g14/taffet2/graph/custom_part.csv"
        }
    };
struct svr_msg
{
    enum svr_event svr_event_type;
    tw_lpid src;          /* source of this request or ack */
    int incremented_flag; /* helper for reverse computation */
    model_net_event_return event_rc;
};

struct {
    int n1; /* The width of the grid */
    int n2; /* The width of each sub-block of the grid */
    int n3; /* The height of each sub-block of the grid */
} cart_info = {0};
int *rand_node_map;
int *rand_node_reverse_map;
static void gen_random_mappings();
static void gen_hardcoded_mappings();
static void rowmaj_lintoxy(int lin_id, int* outx, int* outy) ;
static void blocked_lintoxy(int lin_id, int* outx, int* outy) ;
static void rand_lintoxy(int lin_id, int* outx, int* outy) ;
static int rowmaj_xytolin(int x, int y) ;
static int blocked_xytolin(int x, int y) ;
static int rand_xytolin(int x, int y) ;
 
static void svr_init(
    svr_state * ns,
    tw_lp * lp);
static void svr_event(
    svr_state * ns,
    tw_bf * b,
    svr_msg * m,
    tw_lp * lp);
static void svr_rev_event(
    svr_state * ns,
    tw_bf * b,
    svr_msg * m,
    tw_lp * lp);
static void svr_finalize(
    svr_state * ns,
    tw_lp * lp);

tw_lptype svr_lp = {
    (init_f) svr_init,
    (pre_run_f) NULL,
    (event_f) svr_event,
    (revent_f) svr_rev_event,
    (commit_f) NULL,
    (final_f)  svr_finalize,
    (map_f) codes_mapping,
    sizeof(svr_state),
};

/* setup for the ROSS event tracing
 */
void ft_svr_event_collect(svr_msg *m, tw_lp *lp, char *buffer, int *collect_flag)
{
    (void)lp;
    (void)buffer;
    (void)collect_flag;

    int type = (int) m->svr_event_type;
    memcpy(buffer, &type, sizeof(type));
}

/* can add in any model level data to be collected along with simulation engine data
 * in the ROSS instrumentation.  Will need to update the last field in 
 * ft_svr_model_types[0] for the size of the data to save in each function call
 */
void ft_svr_model_stat_collect(svr_state *s, tw_lp *lp, char *buffer)
{
    (void)s;
    (void)lp;
    (void)buffer;

    return;
}

st_model_types ft_svr_model_types[] = {
    {(ev_trace_f) ft_svr_event_collect,
     sizeof(int),
     (model_stat_f) ft_svr_model_stat_collect,
     0,
     NULL,
     NULL,
     0},
    {NULL, 0, NULL, 0, NULL, NULL, 0}
};

static const st_model_types  *ft_svr_get_model_stat_types(void)
{
    return(&ft_svr_model_types[0]);
}

void ft_svr_register_model_stats()
{
    st_model_type_register("server", ft_svr_get_model_stat_types());
}

const tw_optdef app_opt [] =
{
        TWOPT_GROUP("Model net synthetic traffic " ),
	TWOPT_UINT("traffic", traffic, "UNIFORM RANDOM=1, BISECTION=2, TO_RANK_0=3, PFISTER=4, GRID=5, CSTENCIL=7. Add 1000 to 4, 5 for a different mapping."),
    	TWOPT_UINT("num_messages", num_msgs, "Number of messages to be generated per terminal "),
    	TWOPT_UINT("payload_sz", payload_size, "Size in bytes of each message to generate"),
	TWOPT_STIME("arrival_time", arrival_time, "INTER-ARRIVAL TIME"),
        TWOPT_STIME("load", load, "percentage of terminal link bandiwdth to inject packets"),
        TWOPT_END()
};

const tw_lptype* svr_get_lp_type()
{
            return(&svr_lp);
}

static void svr_add_lp_type()
{
  lp_type_register("server", svr_get_lp_type());
}

static void issue_event(
    svr_state * ns,
    tw_lp * lp,
    bool extra_delay)
{
    (void)ns;
    tw_event *e;
    svr_msg *m;
    tw_stime kickoff_time;

    /* each server sends a dummy event to itself that will kick off the real
     * simulation
     */

    int this_packet_size = 0;
    double this_link_bandwidth = 0.0;

    configuration_get_value_int(&config, "PARAMS", "packet_size", NULL, &this_packet_size);
    if(!this_packet_size) {
        this_packet_size = 0;
        fprintf(stderr, "Packet size not specified, setting to %d\n", this_packet_size);
        exit(0);
    }

    configuration_get_value_double(&config, "PARAMS", "link_bandwidth", NULL, &this_link_bandwidth);
    if(!this_link_bandwidth) {
        this_link_bandwidth = 4.7;
        fprintf(stderr, "Bandwidth of channels not specified, setting to %lf\n", this_link_bandwidth);
    }

    if(arrival_time!=0)
    {
        MEAN_INTERVAL = arrival_time;
    }
    if(load != 0)
    {
        MEAN_INTERVAL = bytes_to_ns(this_packet_size, load*this_link_bandwidth);
    }

    /* skew each kickoff event slightly to help avoid event ties later on */
//    kickoff_time = 1.1 * g_tw_lookahead + tw_rand_exponential(lp->rng, arrival_time);
    kickoff_time = g_tw_lookahead + tw_rand_exponential(lp->rng, MEAN_INTERVAL);
    if (extra_delay)
        kickoff_time += bytes_to_ns(payload_size, load*this_link_bandwidth);


    e = tw_event_new(lp->gid, kickoff_time, lp);
    m = tw_event_data(e);
    m->svr_event_type = KICKOFF;
    tw_event_send(e);
}

static void svr_init(
    svr_state * ns,
    tw_lp * lp)
{
    ns->start_ts = 0.0;

    if (traffic % MAP_SCALE == GRID) {
        int recv_per_round = 4;
        int local_id = codes_mapping_get_lp_relative_id(lp->gid, 0, 0);
        int curx, cury;
        switch(traffic/MAP_SCALE) {
            case REGULAR:
                rowmaj_lintoxy(local_id, &curx, &cury);
                break;
            case TILED:
                blocked_lintoxy(local_id, &curx, &cury);
                break;
            case RANDOM:
            case RCM:
            case METIS:
            case IMPROVEMENT_PARTITION:
                rand_lintoxy(local_id, &curx, &cury);
                break;
        }
        // Unsigned distances to the 4 boundaries
        int distances[] = {curx, cart_info.n1 - 1 - curx, cury, num_nodes/cart_info.n1 - 1 - cury};
        for (int j = 0; j < 4; j++) {
            if (distances[j] == 0)
                recv_per_round--;
            else if (distances[j] == 1)
                recv_per_round++;
        }
        ns->recv_per_round = recv_per_round;
    }

    for (int j = 0; j < 1; j++)
        issue_event(ns, lp, false);
    return;
}

static void handle_kickoff_rev_event(
            svr_state * ns,
            tw_bf * b,
            svr_msg * m,
            tw_lp * lp)
{
    (void)b;
    (void)m;
	ns->msg_sent_count--;
	model_net_event_rc2(lp, &m->event_rc);
    tw_rand_reverse_unif(lp->rng);
}	
static void gen_random_mappings() {
    // Deterministically create a random-looking map
    rand_node_map = malloc(sizeof(int) * num_nodes);
    rand_node_reverse_map = malloc(sizeof(int) * num_nodes);
    int prim_root, prime;
    if (num_nodes == 3564) {
        prim_root = 1275;
        prime = 3571;
    } else if (num_nodes == 4608) {
        prim_root = 2411;
        prime = 4621;
    } else if (num_nodes == 3072) {
        prime = 3079;
        prim_root = 2727;
    }
    int i = 0;
    int acc = 1;
    while (i < num_nodes) {
        acc = (acc * prim_root) % prime;
        if (acc <= num_nodes) {
            rand_node_map[i] = acc - 1;
            rand_node_reverse_map[acc - 1] = i;
            i++;
        }
    }

}
static void gen_hardcoded_mappings() {
    rand_node_map = malloc(sizeof(int) * num_nodes);
    rand_node_reverse_map = malloc(sizeof(int) * num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        if (traffic / MAP_SCALE == RCM) {
            rand_node_map[i] = hardcoded_mappings[traffic / MAP_SCALE][i] - 1;
            rand_node_reverse_map[rand_node_map[i]] = i;
        }
        else {
            rand_node_reverse_map[i] = hardcoded_mappings[traffic / MAP_SCALE][i];
            rand_node_map[rand_node_reverse_map[i]] = i;
        }
    }
}
static void rowmaj_lintoxy(int lin_id, int* outx, int* outy) {
    *outx = lin_id / (num_nodes / cart_info.n1);
    *outy = lin_id % (num_nodes / cart_info.n1);
}
static void blocked_lintoxy(int lin_id, int* outx, int* outy) {
    int lin_blockn = lin_id / (cart_info.n2*cart_info.n3);
    *outx = (lin_blockn % (cart_info.n1/cart_info.n2)) * cart_info.n2;
    *outy = lin_blockn / (cart_info.n1/cart_info.n2) * cart_info.n3;
    int block_off = lin_id % (cart_info.n2*cart_info.n3);
    *outx += block_off / cart_info.n3;
    *outy += block_off % cart_info.n3;
}
static void rand_lintoxy(int lin_id, int* outx, int* outy) {
    rowmaj_lintoxy(rand_node_map[lin_id], outx, outy);
}
static int rowmaj_xytolin(int x, int y) {
    return x * (num_nodes / cart_info.n1) + y;
}
static int blocked_xytolin(int x, int y) {
    int blockx = x/cart_info.n2;
    int blocky = y/cart_info.n3;
    int lin_blockn = blockx + blocky * cart_info.n1/cart_info.n2;
    int lin_block_off = (x % cart_info.n2)*cart_info.n3 + (y % cart_info.n3);
    return lin_blockn * cart_info.n2 * cart_info.n3 + lin_block_off;
}
static int rand_xytolin(int x, int y) {
    return rand_node_reverse_map[rowmaj_xytolin(x,y)];
}
static void handle_kickoff_event(
	    svr_state * ns,
	    tw_bf * b,
	    svr_msg * m,
	    tw_lp * lp)
{
    (void)b;
    (void)m;
//    char* anno;
    char anno[MAX_NAME_LENGTH];
    tw_lpid local_dest = -1, global_dest = -1;

    if(ns->msg_sent_count >= num_msgs)
    {
        m->incremented_flag = 1;
        return;
    }

    m->incremented_flag = 0;
   codes_mapping_get_lp_info(lp->gid, group_name, &group_index, lp_type_name, &lp_type_index, anno, &rep_id, &offset);
   int local_id = codes_mapping_get_lp_relative_id(lp->gid, 0, 0);
   
    svr_msg * m_local = malloc(sizeof(svr_msg));
    svr_msg * m_remote = malloc(sizeof(svr_msg));

    m_local->svr_event_type = LOCAL;
    m_local->src = lp->gid;

    memcpy(m_remote, m_local, sizeof(svr_msg));
    m_remote->svr_event_type = REMOTE;

    ns->start_ts = tw_now(lp);

    switch(traffic % 1000) {
        case UNIFORM:
            {
                /* in case of uniform random traffic, send to a random destination. */
                local_dest = tw_rand_integer(lp->rng, 0, num_nodes - 1);
            }
            break;
        case RANK_ZERO:
            {
                if (local_id != 0)
                    local_dest = 0;
                else {
                    // Rank 0 doesn't send
                    free(m_local);
                    free(m_remote);
                    return;
                }
            }
            break;
        case BISECTION:
            {
                local_dest = (local_id + num_nodes/2) % num_nodes;
                if (traffic / MAP_SCALE == RANDOM) {
                    // Use random permutation to map first half to second half
                    // and so that if x sends to y, then y also sends to x.
                    // The way rand_node_map and rand_node_reverse_map are constructed
                    // means that these loops should terminate
                    if (local_id < num_nodes/2) {
                        local_dest = local_id;
                        do {
                            local_dest = rand_node_map[local_dest];
                        } while (local_dest >= num_nodes / 2);
                        local_dest += num_nodes/2;
                    } else {
                        local_dest = local_id - num_nodes/2;
                        do {
                            local_dest = rand_node_reverse_map[local_dest];
                        } while (local_dest >= num_nodes / 2);
                    }
                }
            }
            break;
        case PFISTER:
            {
                bool in_set; // Part of the 20% of ranks that contribute to congestion?
                if (traffic / 1000 == 0) {
                    local_dest = num_nodes/5 + 1 + tw_rand_integer(lp->rng, 0, 4*num_nodes/5 - 1);
                    in_set = local_id > 0 && local_id <= 712;
                } else {
                    int rand = tw_rand_integer(lp->rng, 0, 4*num_nodes/5 - 1);
                    local_dest = 5*(rand/4);
                    if (rand % 4 > 0)
                        local_dest += 1 + (rand % 4);
                    in_set = local_id % 5 == 1;
                }
                if ((tw_rand_integer(lp->rng, 0, 4) != 0) && in_set)
                {
                    local_dest =  0;
                }
            }
            break;
        case GRID: 
            {
                int curx, cury;
                switch(traffic/MAP_SCALE) {
                    case REGULAR:
                        rowmaj_lintoxy(local_id, &curx, &cury);
                        break;
                    case TILED:
                        blocked_lintoxy(local_id, &curx, &cury);
                        break;
                    case RANDOM:
                    case RCM:
                    case METIS:
                    case IMPROVEMENT_PARTITION:
                        rand_lintoxy(local_id, &curx, &cury);
                        break;
                }
                int destx, desty;
                int direction = tw_rand_integer(lp->rng, 0, 3);
                direction = (ns->msg_sent_count) % 4;
                if (direction < 2) {
                    destx = curx + 2*direction - 1;
                    desty = cury;
                }
                else {
                    desty = cury + 2*direction - 5;
                    destx = curx;
                }

                if (destx < 0 || desty < 0 || destx >= cart_info.n1 || desty >= num_nodes / cart_info.n1) {
                   free(m_local);
                   free(m_remote);
                   ns->msg_sent_count++;
                   printf("Message from %d would be out of bounds. Waiting for next message\n", local_id);
                   issue_event(ns, lp, true);
                   return;
                }
                /*if (destx < 0)
                    destx *= -1;
                if (desty < 0)
                    desty *= -1;
                if (destx >= cart_info.n1)
                    destx = 2*(cart_info.n1 - 1) - destx;
                if (desty >= num_nodes/cart_info.n1)
                    desty = 2*(num_nodes/cart_info.n1 - 1) - desty;*/
                switch(traffic/MAP_SCALE) {
                    case REGULAR:
                        local_dest = rowmaj_xytolin(destx, desty);
                        break;
                    case TILED:
                        local_dest = blocked_xytolin(destx, desty);
                        break;
                    case RANDOM:
                    case RCM:
                    case METIS:
                    case IMPROVEMENT_PARTITION:
                        local_dest = rand_xytolin(destx, desty);
                        break;
                }
            }
            break;
        case UNIFORM_NN:
            {
                int pc = __builtin_popcount(local_id*local_id);
                if (pc % 2 == 0) {
                    if (traffic / MAP_SCALE == 2) {
                        free(m_local);
                        free(m_remote);
                        ns->msg_sent_count++;
                        issue_event(ns, lp, true);
                        return;
                    }
                    if (local_id == 0)
                        local_dest = num_nodes - 1;
                    else
                        local_dest = local_id - 1;
                    while (__builtin_popcount(local_dest*local_dest) % 2 != 0)
                        local_dest--;
                }
                else {
                    do {
                        if (traffic / MAP_SCALE == 1) {
                            free(m_local);
                            free(m_remote);
                            ns->msg_sent_count++;
                            issue_event(ns, lp, true);
                            return;
                        }
                        local_dest = tw_rand_integer(lp->rng, 0, num_nodes - 1);
                    }
                    while(local_dest != local_id && __builtin_popcount(local_dest*local_dest) % 2 == 0);
                }
            }
        case CSTENCIL:
            {
                assert(num_nodes == 3072);
                int dims[3] = {16, 16, 12};
                int direction = (ns->msg_sent_count) % 6;
                int pos[3] = {0, 0, 0};
                int src = local_id;
                for (int i = 0; i < 3; i++) {
                    pos[i] = src % dims[i];
                    src /= dims[i];
                }
                int delta = 2*(direction % 2) - 1;
                pos[direction / 2] = (pos[direction / 2] + delta + dims[direction/2]) % dims[direction/2];

                int scale = 1;
                local_dest = 0;
                for (int i = 0; i < 3; i++) {
                    local_dest += pos[i] * scale;
                    scale *= dims[i];
                }

            }
                

    }
   // printf("Sending from %d to %ld\n", local_id, local_dest);

   assert(local_dest < LLU(num_nodes));

   global_dest = codes_mapping_get_lpid_from_relative(local_dest, group_name, lp_type_name, NULL, 0);

   // printf("global_src:%d, local_src:%d, global_dest:%d, local_dest:%d num_nodes:%d \n",(int)lp->gid, local_id, (int)global_dest,(int)local_dest, num_nodes);

   // If Destination is self, then generate new destination
   if((int)global_dest == (int)lp->gid)
   {
       // printf("Destination was self. Regenerating\n");
       local_dest = (local_dest+1) % (num_nodes-1);
       global_dest = codes_mapping_get_lpid_from_relative(local_dest, group_name, lp_type_name, NULL, 0);
   }

   ns->msg_sent_count++;

   m->event_rc = model_net_event(net_id, "test", global_dest, payload_size, 0.0, sizeof(svr_msg), (const void*)m_remote, sizeof(svr_msg), (const void*)m_local, lp);

   //printf("LP:%d localID:%d Here\n",(int)lp->gid, (int)local_dest);
   //printf("Just Checking net_id:%d\n",net_id);
   // issue_event(ns, lp);
   return;
}

static void handle_remote_rev_event(
            svr_state * ns,     
            tw_bf * b,
            svr_msg * m,
            tw_lp * lp)
{
        (void)b;
        (void)m;
        (void)lp;
        ns->msg_recvd_count--;
}

static void handle_remote_event(
	    svr_state * ns,
	    tw_bf * b,
	    svr_msg * m,
	    tw_lp * lp)
{
    (void)b;
    (void)m;
    (void)lp;
    ns->msg_recvd_count++;

    if (traffic % 1000 == GRID) {
        if (ns->msg_recvd_count % ns->recv_per_round == 0) {
        printf("Node %d finished round %d at time %f\n", codes_mapping_get_lp_relative_id(lp->gid, 0, 0), ns->msg_recvd_count / ns->recv_per_round, tw_now(lp));
            for (int i = 0; i < 0; i++)
                issue_event(ns, lp, false);
        }
    } 
    else if (traffic % 1000 == CSTENCIL)
    {
    }
    else {
        tw_event *e;
        e = tw_event_new(m->src, .01, lp);
        svr_msg * ma = tw_event_data(e);
        ma->svr_event_type = ACK;
        tw_event_send(e);
    }
}
static void handle_ack_rev_event(
	    svr_state * ns,
	    tw_bf * b,
	    svr_msg * m,
	    tw_lp * lp)
{
    (void)b;
    (void)m;
    (void)lp;
    printf("How do I handle this reverse event??!\n");
    assert(0);
}
static void handle_ack_event(
	    svr_state * ns,
	    tw_bf * b,
	    svr_msg * m,
	    tw_lp * lp)
{
    (void)b;
    (void)m;
    (void)lp;
   issue_event(ns, lp, false);
}

static void handle_local_rev_event(
                svr_state * ns,
                tw_bf * b,
                svr_msg * m,
                tw_lp * lp)
{
    (void)b;
    (void)m;
    (void)lp;
	ns->local_recvd_count--;
}

static void handle_local_event(
                svr_state * ns,
                tw_bf * b,
                svr_msg * m,
                tw_lp * lp)
{
    (void)b;
    (void)m;
    (void)lp;
    ns->local_recvd_count++;
    if (traffic % 1000 == GRID || traffic % 1000 == CSTENCIL) {
        issue_event(ns, lp, false);
    }
}

static void svr_finalize(
    svr_state * ns,
    tw_lp * lp)
{
    ns->end_ts = tw_now(lp);

//    printf("server %llu recvd %d bytes in %f seconds, %f MiB/s sent_count %d recvd_count %d local_count %d \n", (unsigned long long)lp->gid, PAYLOAD_SZ*ns->msg_recvd_count, ns_to_s(ns->end_ts-ns->start_ts),
//        ((double)(PAYLOAD_SZ*ns->msg_sent_count)/(double)(1024*1024)/ns_to_s(ns->end_ts-ns->start_ts)), ns->msg_sent_count, ns->msg_recvd_count, ns->local_recvd_count);
    return;
}

static void svr_rev_event(
    svr_state * ns,
    tw_bf * b,
    svr_msg * m,
    tw_lp * lp)
{
    switch (m->svr_event_type)
    {
	case REMOTE:
		handle_remote_rev_event(ns, b, m, lp);
		break;
	case LOCAL:
		handle_local_rev_event(ns, b, m, lp);
		break;
	case KICKOFF:
		handle_kickoff_rev_event(ns, b, m, lp);
		break;
	case ACK:
		handle_ack_rev_event(ns, b, m, lp);
		break;
	default:
		assert(0);
		break;
    }
}

static void svr_event(
    svr_state * ns,
    tw_bf * b,
    svr_msg * m,
    tw_lp * lp)
{
   switch (m->svr_event_type)
    {
        case REMOTE:
            handle_remote_event(ns, b, m, lp);
            break;
        case LOCAL:
            handle_local_event(ns, b, m, lp);
            break;
	case KICKOFF:
	    handle_kickoff_event(ns, b, m, lp);
	    break;
        case ACK:
            handle_ack_event(ns, b, m, lp);
            break;
        default:
            printf("\n LP: %d has received invalid message from src lpID: %d of message type:%d", (int)lp->gid, (int)m->src, m->svr_event_type);
            assert(0);
        break;
    }
}

int main(
    int argc,
    char **argv)
{

    int nprocs;
    int rank;
    int num_nets;
    int *net_ids;

    lp_io_handle handle;

    tw_opt_add(app_opt);

    tw_init(&argc, &argv);
#ifdef USE_RDAMARIS
    if(g_st_ross_rank)
    { // keep damaris ranks from running code between here up until tw_end()
#endif
    codes_comm_update();

    offset = 1;

    if(argc < 2)
    {
            printf("\n Usage: mpirun <args> --sync=2/3 mapping_file_name.conf (optional --nkp) ");
            MPI_Finalize();
            return 0;
    }

    MPI_Comm_rank(MPI_COMM_CODES, &rank);
    MPI_Comm_size(MPI_COMM_CODES, &nprocs);

    configuration_load(argv[2], MPI_COMM_CODES, &config);

    model_net_register();

    svr_add_lp_type();

    if (g_st_ev_trace || g_st_model_stats || g_st_use_analysis_lps)
        ft_svr_register_model_stats();

    codes_mapping_setup();


    net_ids = model_net_configure(&num_nets);
    //assert(num_nets==1);
    net_id = *net_ids;
    free(net_ids);

    if(net_id != FATTREE)
    {
	printf("\n The test works with fat tree model configuration only! ");
        MPI_Finalize();
        return 0;
    }
    num_servers_per_rep = codes_mapping_get_lp_count("MODELNET_GRP", 1, "server",
            NULL, 1);
    configuration_get_value_int(&config, "PARAMS", "num_routers", NULL, &num_routers_per_grp);
    
    num_groups = (num_routers_per_grp * (num_routers_per_grp/2) + 1);
    num_nodes = num_groups * num_routers_per_grp * (num_routers_per_grp / 2);
    num_nodes_per_grp = num_routers_per_grp * (num_routers_per_grp / 2);

    num_nodes = codes_mapping_get_lp_count("MODELNET_GRP", 0, "server", NULL, 1);

    printf("num_nodes:%d \n",num_nodes);

    cart_info.n1 = (int) sqrt(num_nodes);
    while (num_nodes % cart_info.n1 > 0)
        cart_info.n1++;
    // TODO: Don't hardcode this
    if (num_nodes == 3564) {
        cart_info.n2 = 3;
        cart_info.n3 = 6;
    }
    if (num_nodes == 4608) {
        cart_info.n2 = 8;
        cart_info.n3 = 4;
    }

    if (traffic / 1000 == RANDOM) {
        gen_random_mappings();
    }
    if (traffic / 1000 >= RCM) {
        gen_hardcoded_mappings();
    }


    if(lp_io_prepare("modelnet-test", LP_IO_UNIQ_SUFFIX, &handle, MPI_COMM_WORLD) < 0)
    {
        return(-1);
    }
    modelnet_stats_dir = lp_io_handle_to_dir(handle);

    tw_run();

    model_net_report_stats(net_id);


#if PARAMS_LOG
    if(!g_tw_mynode)
    {
	char temp_filename[1024];
	char temp_filename_header[1024];
	sprintf(temp_filename,"%s/sim_log.txt",modelnet_stats_dir);
	sprintf(temp_filename_header,"%s/sim_log_header.txt",modelnet_stats_dir);
	FILE *fattree_results_log=fopen(temp_filename, "a");
	FILE *fattree_results_log_header=fopen(temp_filename_header, "a");
	if(fattree_results_log == NULL)
		printf("\n Failed to open results log file %s in synthetic-fattree\n",temp_filename);
	if(fattree_results_log_header == NULL)
		printf("\n Failed to open results log header file %s in synthetic-fattree\n",temp_filename_header);
	printf("Printing Simulation Parameters/Results Log File\n");
	fprintf(fattree_results_log_header,", <Workload>, <Load>, <Mean Interval>, ");
	fprintf(fattree_results_log,"%11.3d, %5.2f, %15.2f, ",traffic, load, MEAN_INTERVAL);
	fclose(fattree_results_log_header);
	fclose(fattree_results_log);
    }
#endif

    if(lp_io_flush(handle, MPI_COMM_CODES) < 0)
    {
        return(-1);
    }
#ifdef USE_RDAMARIS
    } // end if(g_st_ross_rank)
#endif
    tw_end();

#if PARAMS_LOG
    if(!g_tw_mynode)
    {
	char temp_filename[1024];
	char temp_filename_header[1024];
	sprintf(temp_filename,"%s/sim_log.txt",modelnet_stats_dir);
	sprintf(temp_filename_header,"%s/sim_log_header.txt",modelnet_stats_dir);
	FILE *fattree_results_log=fopen(temp_filename, "a");
	FILE *fattree_results_log_header=fopen(temp_filename_header, "a");
	FILE *fattree_ross_csv_log=fopen("ross.csv", "r");
	if(fattree_results_log == NULL)
		printf("\n Failed to open results log file %s in synthetic-fattree\n",temp_filename);
	if(fattree_results_log_header == NULL)
		printf("\n Failed to open results log header file %s in synthetic-fattree\n",temp_filename_header);
	if(fattree_ross_csv_log == NULL)
		tw_error(TW_LOC, "\n Failed to open ross.csv log file \n");
	printf("Reading ROSS specific data from ross.csv and Printing to Fat Tree Log File\n");
	
	char * line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, fattree_ross_csv_log);
	while (read != -1) 
	{
		read = getline(&line, &len, fattree_ross_csv_log);
	}

	char * pch;
	pch = strtok (line,",");
	int idx = 0;
        int gvt_computations;
	long long total_events, rollbacks, net_events;
        float running_time, efficiency, event_rate;
	while (pch != NULL)
	{
		pch = strtok (NULL, ",");
		switch(idx)
		{
			case 4:
				total_events = atoll(pch);
				break;
			case 13:
				rollbacks = atoll(pch);
				break;
			case 17:
				gvt_computations = atoi(pch);
				break;
			case 18:
				net_events = atoll(pch);
				break;
			case 3:
				running_time = atof(pch);
				break;
			case 8:
				efficiency = atof(pch);
				break;
			case 19:
				event_rate = atof(pch);
				break;
		}
		idx++;
	}
	fprintf(fattree_results_log_header,"<Total Events>, <Rollbacks>, <GVT Computations>, <Net Events>, <Running Time>, <Efficiency>, <Event Rate>");
	fprintf(fattree_results_log,"%14llu, %11llu, %18d, %12llu, %14.4f, %12.2f, %12.2f\n",total_events,rollbacks,gvt_computations,net_events,running_time,efficiency,event_rate);
	fclose(fattree_results_log);
	fclose(fattree_results_log_header);
	fclose(fattree_ross_csv_log);
    }
#endif

    return 0;
}

/*
 * Local variables:
 *  c-indent-level: 4
 *  c-basic-offset: 4
 * End:
 *
 * vim: ft=c ts=8 sts=4 sw=4 expandtab
 */
