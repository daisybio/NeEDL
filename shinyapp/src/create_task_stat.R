###############################
#### Create Task Statistic ####
###############################


#' Create task statistics (from LRZ file)
#'
#' @param task_stat_path Path to task statistics file
#'
#' @return Parsed task statistics file 
#'
#' @examples create_task_stats(tast_stat_path)
create_task_stats <- function(task_stat_path){
  task_stat_list <- readLines(task_stat_path)
  
  # Command 
  command_dt <- data.table(unlist(str_split(task_stat_list[1], "--")))
  command_dt$V1 <- gsub("\tCommand being timed:", "Command_being_timed", command_dt$V1)
  command_dt$V1 <- gsub('"', '', command_dt$V1)
  command_dt <- separate(command_dt, col = "V1", into = c("active flag", "value"), sep = " ")
  
  # Other stats in file 
  other_stats_dt <- data.table(task_stat_list[-1])
  other_stats_dt <- separate(other_stats_dt, col = "V1", into = c("active flag", "value"), sep = ": ")
  
  # Merged tast stats data table
  task_stat_dt <- rbind(command_dt, other_stats_dt)
  
  return(task_stat_dt)
}


#' Create statistics for time and memory
#'
#' @param run_log_path Path to run.log file 
#'
#' @return Parsed data table with time and memory information
#'
#' @examples create_stats_from_run_log(run_log_path)
create_stats_from_run_log <- function(run_log_path){
  
  # read run.log file
  run_log_list <- readLines(run_log_path)
  
  
  # if operating system is not Linux the RUN STATS are not in the run.log file
  if(sum(grepl("RUN STATS:", run_log_list)) == 0){
    run_log_dt <- data.table(statistic = c("Not available"), value = c("Not available"))
    return(run_log_dt)
  }
  
  # get index where run stats start + 1
  run_stats_start <- match(TRUE, grepl("RUN STATS:", run_log_list)) + 1
  
  # get entries starting from one line after run stats start
  run_log_list <- run_log_list[run_stats_start:length(run_log_list)]
  
  # remove time information 
  run_log_list <- unlist(lapply(run_log_list, function(x) unlist(str_split(x, "]      "))[2]))
  
  # Create data table with statistic - value information
  run_log_dt <- separate(data.table(run_log_list), col = "run_log_list", into = c("statistic", "value"), sep = ": ")
  
  
  return(run_log_dt)
}





#' Create parameter overview from pipeline config json file
#'
#' @param pipeline_config_json Pipeline config json file
#'
#' @return Data table containing all jobs information from the json (not JOB: sequential)
#'
#' @examples create_param_from_pipeline_config(pipeline_config_json)
create_param_from_pipeline_config <- function(pipeline_config_json){
  # Parse job 1 to 4
  jobs_json <- data.table(pipeline_config_json$jobs %>% spread_all())
  
  jobs_1_to_4_dt <- rbindlist(lapply(seq(1, dim(jobs_json)[1] - 1), function(job_index){
    tmp_job <- (unlist(jobs_json[job_index]$..JSON))
    tmp_job_dt <- data.table(parameter = names(tmp_job), values = tmp_job)
    tmp_job_dt
  }))
  
  
  # Parse job 5
  job_5_json <- data.table(pipeline_config_json$jobs[[5]]$jobs[[1]]$config$jobs %>% spread_all())
  
  job_5_dt <- rbindlist(lapply(seq(1, dim(job_5_json)[1]), function(job_index){
    tmp_job <- (unlist(job_5_json[job_index]$..JSON))
    tmp_job_dt <- data.table(parameter = names(tmp_job), values = tmp_job)
    tmp_job_dt
  }))
  
  
  # Merge job 1 to 4 with 5
  all_jobs_dt <- rbind(jobs_1_to_4_dt, job_5_dt)
  
  num_jobs <- dim(all_jobs_dt[parameter == "JOB"])[1]
  
  all_jobs_dt[parameter == "JOB"]$parameter <- paste0("<strong>","Pipeline step ", seq(1, num_jobs), "/", num_jobs, "</strong>")
  
  
  return(all_jobs_dt)
}





