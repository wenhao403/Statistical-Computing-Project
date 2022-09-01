library(shiny)
library(dplyr)
library(purrr)
library(ggplot2)
library(gganimate)
library(plotly)
library(shinythemes)

library(coda)
library(forecast)
library(patchwork)

## 2022.9.1
library(gifski)
library(rsconnect)


########  UI  ########
ui <- fluidPage(
  navbarPage('Statistical Computing Final Project',
             
             ########  tabPanel 1  ########
             tabPanel('Density Plot', fluidPage(theme=shinytheme('flatly')),
                      
                      ##  sidebarLayout
                      sidebarLayout(
                        sidebarPanel(width=3,
                          sliderInput('mu', 'Mu :', min=-10, max=10, step=0.1, value=2),
                          sliderInput('sigma', 'Sigma :', min=1, max=7, step=0.1, value=3),
                          # h5('plot density : '),
                          actionButton('go1_1', 'Go')),
                        
                        ##  mainPanel
                        mainPanel(
                          tabsetPanel(
                            ########  Tab 1  ########
                            tabPanel('Density function',
                                     h3(),
                                     withMathJax(),
                                     splitLayout(cellWidths=500, 
                                                 plotOutput("true_curve1_1", height=400),
                                                 plotOutput("true_curve1_2", height=400)),
                                     uiOutput('ex1_1')),
                            
                            ########  Tab 2  ########
                            tabPanel('About',
                                     h3(),
                                     uiOutput('ex1_2')))
                          )
                        )
                      ),
             
             ########  tabPanel 2  ########
             tabPanel('Metropolis Hastings',
                      
                      ##  sidebarLayout
                      sidebarLayout(
                        sidebarPanel(width=3,
                                     numericInput('n', 'MH sample size (a) :', value=20000),
                                     actionButton('go2_1', 'Go1'),
                                     h3(),
                                     numericInput('gelman_plot_intercept_number', 'Gelman plot converage cut off (b) :', value=1000),
                                     numericInput('acf_max_lag', 'ACF max lag :', value=100),
                                     numericInput('acf_lag', 'ACF lag cut off (c) :', value=5),
                                     actionButton('go2_2', 'Go2'),
                                     h3(),
                                     radioButtons("static_dynamic", "Static or Dynamic",
                                                  choices = c('Static' = 'plot',
                                                              'Dynamic' = 'plot2'),
                                                  selected = 'plot'),
                                     h5('Help : Final sample size = ( a - 4b ) / c')),
                        
                        ##  mainPanel
                        mainPanel(
                          tabsetPanel(
                            ########  Tab 1  ########
                            tabPanel('Gelman & ACF plot',
                                     h4('Gelman plot'),
                                     plotOutput('gelman_plot'),
                                     h4('ACF plot'),
                                     plotOutput('acf1')),
                            
                            ########  Tab 2  ########
                            tabPanel('Sampling data',
                                     conditionalPanel(
                                       condition = "input.static_dynamic == 'plot'",
                                       h3(),
                                       splitLayout(cellWidths=700, plotOutput("plot", height=400)),
                                       tableOutput('values2_1'),
                                       tableOutput('values2_2')),
                                     conditionalPanel(
                                       condition = "input.static_dynamic == 'plot2'",
                                       h3(),
                                       imageOutput('plot2'),
                                       tableOutput('values_dynamic_hist2_1'),
                                       tableOutput('values_dynamic_hist2_2'))),
                            
                            ########  Tab 3  ########
                            tabPanel('About',
                                     h3(),
                                     uiOutput('ex2_1')))
                          )
                        )
                      ),
             
             ########  tabPanel 3  ########
             tabPanel('Coordiante Descent',
                      
                      ##  sidebarLayout
                      sidebarLayout(
                        sidebarPanel(width=3,
                                     sliderInput('init_mu', 'Initial Mu :', min=-10, max=10, step=0.1, value=1),
                                     sliderInput('init_sigma', 'Initial Sigma :', min=1, max=5, step=0.1, value=2.5),
                                     numericInput('n_iter', 'Iterative Number :', value=100),
                                     # h5('Interative :')
                                     actionButton('go3_1', 'Go'),
                                     h3(),
                                     radioButtons("static_dynamic2", "Static or Dynamic",
                                                  choices = c('Static' = 'plot_cd_para_static',
                                                              'Dynamic' = 'plot_cd_para'),
                                                  selected = 'plot_cd_para_static')),
                        
                        ##  mainPanel
                        mainPanel(
                          tabsetPanel(
                            ########  Tab 1  ########
                            tabPanel('Parameter estimate',
                                     conditionalPanel(
                                       condition = "input.static_dynamic2 == 'plot_cd_para_static'",
                                       h3(),
                                       splitLayout(cellWidths=500, 
                                                   plotOutput("plot_cd_para_static", height=500),
                                                   plotOutput("plot_cd_loss_static", height=500)),
                                       tableOutput('values3_1')),
                                     conditionalPanel(
                                       condition = "input.static_dynamic2 == 'plot_cd_para'",
                                       h3(),
                                       splitLayout(cellWidths=500,
                                                   imageOutput("plot_cd_para", height=500),
                                                   imageOutput("plot_cd_loss", height=500)),
                                       tableOutput('values3_2'))),
                            
                            ########  Tab 2  ########
                            tabPanel('Log-likelihood',
                                     h3('Log-Likelihood'),
                                     plotlyOutput('plot_ll_3d'),
                                     uiOutput('ex3_1'),
                                     plotlyOutput('plot_ll_heatmap')),
                            
                            ########  Tab 3  ########
                            tabPanel('About',
                                     h3(),
                                     helpText('Log-likelihood function of \\( (\\mu,\\sigma) \\): '),
                                     uiOutput('ex3_2')))
                          )
                        )
                      )
             )
  )






##  global function  
f = function(x, mu, sigma){
  x[x>=0] = 1/(sigma*sqrt(2*pi)) * (exp((-1/2)*((x[x>=0]+mu)/sigma)^2) + exp((-1/2)*((x[x>=0]-mu)/sigma)^2))
  x[x<0] = 0
  return(x)
}
g = function(x, a, b){
  1 / (b-a)  # U(a,b)
}
true_mean = function(mu, sigma){
  sigma * sqrt(2/pi) * exp(-mu^2/(2*sigma^2)) + mu * (1-2*pnorm(-mu/sigma))
}
true_var = function(mu, sigma){
  mu^2 + sigma^2 - (true_mean(mu, sigma))^2
}
true_sd = function(mu, sigma){
  sqrt(mu^2 + sigma^2 - (true_mean(mu, sigma))^2)
}
ll = function(para, input){
  mu = para[1]
  sigma = para[2]
  n = length(input)
  # negative loglikelihood
  - ((-n/2) * log(2*pi*sigma^2) + sum(log(exp((-1/(2*sigma^2))*(input-mu)^2) + exp((-1/(2*sigma^2))*(input+mu)^2)))) 
}

Sampling_MH = function(init,mu,sigma, n){
  x = c()
  succ = c()
  x[1] = init
  set.seed(init)
  for(i in 2:n) {
    x[i] = rlnorm(1,abs(mu),sigma) 
    r = (f(x[i],mu,sigma) / f(x[i-1],mu,sigma)) * (dlnorm(x[i-1],abs(mu),sigma) / dlnorm(x[i],abs(mu),sigma))
    succ[i] <- rbinom(n=1, size=1, prob=min(r,1)) #Ber(p=min(r,1))
    if(succ[i]==0) x[i] = x[i-1] #not move
  }
  return(x)
}


########  SERVER  ########
server <- function(input, output){
    
    ##  variables  
    # react_ = reactiveValues(n=input$n)
    # observe({
    #   mc1 <<- mcmc(Sampling_MH(1, input$mu, input$sigma, (input$n)/4))
    #   mc2 <<- mcmc(Sampling_MH(2, input$mu, input$sigma, (input$n)/4))
    #   mc3 <<- mcmc(Sampling_MH(3, input$mu, input$sigma, (input$n)/4))
    #   mc4 <<- mcmc(Sampling_MH(4, input$mu, input$sigma, (input$n)/4))
    # })
  
    react = reactive({
      # MH
      mc1 = mcmc(Sampling_MH(1, input$mu, input$sigma, input$n/4))
      mc2 = mcmc(Sampling_MH(2, input$mu, input$sigma, input$n/4))
      mc3 = mcmc(Sampling_MH(3, input$mu, input$sigma, input$n/4))
      mc4 = mcmc(Sampling_MH(4, input$mu, input$sigma, input$n/4))
      mc1234 = mcmc.list(mc1, mc2, mc3, mc4)
      
      tmp_x = c(mc1[-(1:input$gelman_plot_intercept_number)],
                mc2[-(1:input$gelman_plot_intercept_number)], 
                mc3[-(1:input$gelman_plot_intercept_number)],
                mc4[-(1:input$gelman_plot_intercept_number)]) # 預設 2000
      
      x = tmp_x[seq(1, length(tmp_x), by=input$acf_lag)]
      
      df_fold_norm_rs = data.frame(sample=1:length(x), x=x)
      
      df_ani = df_fold_norm_rs %>%
        split(.$sample) %>%
        accumulate(~ bind_rows(.x, .y)) %>%
        bind_rows(.id='frame') %>%
        mutate(frame=as.integer(frame))
      
      # h: histogram
      h = hist(df_fold_norm_rs$x, 30, plot=F)
      
      true_mean = true_mean(input$mu, input$sigma)
      est = mean(x)
      est_se = sqrt((1/length(x)) * (1/(length(x)-1))*sum((x-mean(x))^2)) # s.e.(estimator)
      
      true_var = true_var(input$mu, input$sigma)
      y = (x - mean(x))^2
      est_var = mean(y)
      est_var_se = sqrt((1/length(y)) * (1/(length(y)-1))*sum((y-mean(y))^2)) # s.e.(estimator)
      
      ll_fix_sigma = function(input, mu){
        n = length(input$x)
        sigma = input$sigma
        - ((-n/2) * log(2*pi*sigma^2) + sum(log(exp((-1/(2*sigma^2))*(input$x-mu)^2) + exp((-1/(2*sigma^2))*(input$x+mu)^2)))) 
      }
      
      ll_fix_mu = function(input, sigma){
        n = length(input$x)
        mu = input$mu
        - ((-n/2) * log(2*pi*sigma^2) + sum(log(exp((-1/(2*sigma^2))*(input$x-mu)^2) + exp((-1/(2*sigma^2))*(input$x+mu)^2))))
      }
      
      ## coordinte descent
      coordinate.descent = function(init, n.iter){
        para = matrix(NA, n.iter, 2)  # para=(mu, sigma)
        para[1,] = init
        loss = c(ll(para[1,], input=x))
        
        num = c()
        for (i in 2:n.iter){
          # find mu mle
          input = list(x=x, sigma=para[i-1,2])
          mu_hat = nlminb(start=para[i-1,1], objective=ll_fix_sigma, input=input)$par
          input$mu = mu_hat
          
          # find sigma mle
          sigma_hat = nlminb(start=para[i-1,2], objective=ll_fix_mu, input=input)$par
          input$sigma = sigma_hat
          
          para[i,] = c(input$mu, input$sigma)
          loss[i] = c(ll(para[i,], input=input$x))
          
          if (max(abs(para[i,]-para[i-1,])) < 1e-10){
            num = c(num, i)
          }
        }
        return(list(para=para, loss=loss, num=num))
      }
      
      cd = coordinate.descent(init=c(input$init_mu, input$init_sigma), n.iter=input$n_iter)
      cd_para = cd$para
      cd_loss = cd$loss
      
      df_cd_para = data.frame(time=rep(1:input$n_iter, 2), para=c(cd_para[,1], cd_para[,2]),
                           name=c(rep('mu', input$n_iter), rep('sigma', input$n_iter)))
      df_cd_loss = data.frame(time=seq(1,input$n_iter,by=1), loss=cd_loss, name=rep('loss',input$n_iter))
      
      ll_outer = function(mu, sigma){
        n = length(x)
        # loglikelihood
        ((-n/2) * log(2*pi*sigma^2) + sum(log(exp((-1/(2*sigma^2))*(x-mu)^2) + exp((-1/(2*sigma^2))*(x+mu)^2))))
      }
      tmp = c()
      test_x = seq(0, abs(input$mu*2), len=100) # 20: upper limit of mu
      test_y = seq(1, input$sigma*2, len=100) # 5: upper limit of sigma
      for (mu_tmp in test_x){
        for (sigma_tmp in test_y){
          tmp = c(tmp, ll_outer(mu_tmp, sigma_tmp))
        }
      }
      tmp_m = matrix(tmp, length(test_x), length(test_y), byrow=T)
      
      tmp_sample = rnorm(2000, input$mu, input$sigma)
      normal_hist = hist(tmp_sample, plot=F)
      folded_normal_hist = hist(abs(tmp_sample), plot=F)
      density_max = max(normal_hist$density, folded_normal_hist$density)+0.03
      
      react = list(x=tmp_x, df_fold_norm_rs=df_fold_norm_rs, df_ani=df_ani, h=h,
                   true_mean=true_mean, est=est, est_se=est_se,
                   true_var=true_var, est_var=est_var, est_var_se=est_var_se,
                   cd_para=cd_para, cd_loss=cd_loss, df_cd_para=df_cd_para, df_cd_loss=df_cd_loss,
                   test_x=test_x, test_y=test_y, tmp_m=tmp_m,
                   density_max=density_max,
                   mc1234=mc1234, acf_x=x)
      })

    
    ########  tabPanel 1  ########
    ########  Tab1  ########
    output$true_curve1_1 = renderPlot({
      input$go1_1 # values: 0, 1, 2, 3,...
      isolate(
        ggplot(NULL, aes(-40,40)) +
          geom_area(stat='function', fun=dnorm, args=list(input$mu, input$sigma), fill='#F8766D', xlim=c(-40, 0), alpha=0.9) +
          geom_area(stat='function', fun=dnorm, args=list(input$mu, input$sigma), fill='#619CFF', xlim=c(0, 40)) +
          geom_vline(xintercept=0, colour="#F8766D", size=1.2, alpha=0.5) +
          theme_minimal() +
          ggtitle('normal distribution') +
          labs(x='x', y='Density') +
          theme(plot.title=element_text(size=20),
                axis.title=element_text(size=15)) +
          xlim(input$mu-3.5*input$sigma, input$mu+3.5*input$sigma) +
          ylim(c(0, react()$density_max)))
      })
    
    output$true_curve1_2 = renderPlot({
      input$go1_1
      isolate(
        ggplot(NULL, aes(0, 40)) + 
          geom_area(stat="function", fun=f, args=list(abs(input$mu), input$sigma), fill="#F8766D", alpha=0.9) +
          geom_area(stat="function", fun=dnorm, args=list(input$mu, input$sigma), fill="#619CFF") +
          theme_minimal() +
          ggtitle('normal & folded normal distribution') +
          labs(x='x', y='Density') +
          theme(plot.title=element_text(size=20),
                axis.title=element_text(size=15)) +
          xlim(c(0, abs(input$mu)+3.5*input$sigma)) +
          ylim(c(0, react()$density_max)))
      })
    
    output$ex1_1 = renderUI({
      withMathJax(helpText('Folded normal distribution:  $$f(x;\\mu,\\sigma^2)=\\frac{1}{\\sqrt{2\\pi\\sigma^2}}e^{-(x-\\mu)^2/(2\\sigma^2)}+\\frac{1}{\\sqrt{2\\pi\\sigma^2}}e^{-(x+\\mu)^2/(2\\sigma^2)},\ 0<x<\\infty$$'))
      })
    
    ########  Tab2  ########
    output$ex1_2 = renderUI({
      withMathJax(
        helpText('Defintion : '),
        helpText('\\( \\text{Give a normally distributed random variable X with mean } \\mu \\text{ and varianve } \\sigma \\text{, the random variable Y=|X| has a folded normal distribution.}\\)'),
        helpText('Density :'),
        helpText('$$f(x;\\mu,\\sigma^2)=\\frac{1}{\\sqrt{2\\pi\\sigma^2}}e^{-(x-\\mu)^2/(2\\sigma^2)}+\\frac{1}{\\sqrt{2\\pi\\sigma^2}}e^{-(x+\\mu)^2/(2\\sigma^2)},\ 0<x<\\infty,\ -\\infty<\\mu<\\infty,\\sigma>0$$'),
        helpText('CDF :'),
        helpText('$$F(x;\\mu,\\sigma^2)=\\frac{1}{2} (\\text{erf}(\\frac{x+\\mu}{\\sqrt{2\\sigma^2}})+\\text{erf}(\\frac{x-\\mu}{\\sqrt{2\\sigma^2}}))$$'),
        helpText('Mean of folded distribution :'),
        helpText('$$E(X)=\\sigma\\sqrt{\\frac{2}{\\pi}}\\text{exp}(\\frac{-\\mu^{2}}{2\\sigma^{2}})+\\mu\\text{erf}(\\frac{\\mu}{\\sqrt{2\\sigma^2}})$$'),
        helpText('Variance :'),
        helpText('$$Var(X)=\\mu^2+\\sigma^2-E^2(X)$$'))
      })
    
    
    ########  tabPanel 2  ########
    ########  Tab1  ########
    output$gelman_plot = renderPlot({
      input$go1_1
      input$go2_1
      isolate(gelman.plot(react()$mc1234))
      })
    output$acf1 = renderPlot({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate({
        ggAcf(react()$x, lag.max=input$acf_max_lag) +
          ggtitle('ACF plot') +
          ggAcf(react()$acf_x, lag.max=input$acf_max_lag) +
          ggtitle('Lag cut off ACF plot')
        })
      })
    
    ########  Tab2  ########
    ## Static
    output$plot = renderPlot({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate({
        p_anim = react()$df_ani %>%
          ggplot(aes(x=x)) +
          geom_histogram(aes(y=..density..), bins=25, fill="darkgreen", color=1, alpha=0.5) + #fill="#69b3a2"
          stat_function(fun=f, args=list(mu=input$mu, sigma=input$sigma), col='darkgreen', size=2, alpha=1) +
          theme_minimal() +
          theme(plot.title=element_text(size=25),
                axis.title=element_text(size=15),
                legend.title=element_text(size=15),
                legend.text=element_text(size=15)) +
          labs(title='Histogram of data') +
          coord_cartesian(xlim=c(react()$h$breaks[1], react()$h$breaks[length(react()$h$breaks)]),
                          ylim=c(0, max(react()$h$density)*1.2))
        p_anim
        })
      })
    
    ## Dynamic
    output$plot2 = renderImage({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate({
        p_anim = react()$df_ani %>%
          ggplot(aes(x=x)) +
          geom_histogram(aes(y=..density..), bins=25, fill="darkgreen", color=1, alpha=0.5) + #fill="#69b3a2"
          stat_function(fun=f, args=list(mu=input$mu, sigma=input$sigma), col='darkgreen', size=2, alpha=1) +
          # geom_histogram(aes(y=..density..), bins=25, fill="#E7B800", color=1, alpha=0.5) + #fill="#69b3a2"
          # stat_function(fun=f, args=list(mu=input$mu, sigma=input$sigma), col='orange', size=2, alpha=0.8) +
          theme_minimal() +
          theme(plot.title=element_text(size=25),
                axis.title=element_text(size=15),
                legend.title=element_text(size=15),
                legend.text=element_text(size=15)) +
          labs(title='sample:{frame}') +
          coord_cartesian(xlim=c(react()$h$breaks[1], react()$h$breaks[length(react()$h$breaks)]),
                          ylim=c(0, max(react()$h$density)*1.2))
        
        # animate
        anim = p_anim +
          transition_manual(frame)
        anim_save("outfile.gif", animate(anim, renderer=gifski_renderer(loop=T), width=700, height=400)) # New
        # Return a list containing the filename
        list(src="outfile.gif", contentType="image/gif") #, width=1000, height=500
        })
      }, deleteFile = T)
    
    
    sliderValues2_1 = reactive({
      data.frame(
        'Mu' = input$mu,
        'Sigma' = input$sigma,
        'Mean' = react()$true_mean,
        'Var' = react()$true_var,
        'Sample_Size' = length(react()$acf_x))  # 檢查完收斂後的剩餘樣本數
      })
    output$values2_1 <- renderTable({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate(sliderValues2_1())
      })
    output$values_dynamic_hist2_1 = renderTable({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate(sliderValues2_1())
      })
    
    sliderValues2_2 = reactive({
      df2 = data.frame(
        'Est_mean' = c(react()$est, react()$est_se),
        'Est_var' = c(react()$est_var, react()$est_var_se))
      })
    output$values2_2 <- renderTable({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate(sliderValues2_2())
      })
    output$values_dynamic_hist2_2 = renderTable({
      input$go1_1
      input$go2_1
      input$go2_2
      isolate(sliderValues2_2())
      })
    
    ########  Tab3  ########
    output$ex2_1 = renderUI({
      withMathJax(
        helpText('target pdf : $$X \\sim f(x;\\mu,\\sigma^2)=\\frac{1}{\\sqrt{2\\pi\\sigma^2}}e^{-(x-\\mu)^2/(2\\sigma^2)}+\\frac{1}{\\sqrt{2\\pi\\sigma^2}}e^{-(x+\\mu)^2/(2\\sigma^2)},\\ 0<x<\\infty$$'),
        helpText('jump proposal pdf : $$q(x;\\mu,\\sigma^2)=\\frac{1}{x\\sigma\\sqrt{2\\pi}}e^{-(\\text{ln}x-\\mu)^2/(2\\sigma^2)},\\ 0<x<\\infty,\\ \\text{ i.e. } q \\sim \\text{logN}(\\mu, \\sigma^2)$$'),
        # helpText('$$i.e. q \\sim N(|\\mu|, \\sigma^2)$$'),
        helpText('Algorithm : '),
        helpText('\\( \\text{Step1 : Set initial value } x^{(0)}, \\text{ satisfying  } f(x^{(0)})>0 \\)'),
        helpText('\\( \\text{Step2 : For } t=1,2,...., \\)'),
        helpText('\\( \\ \\ \\ \\text{ 1. draw } x^* \\sim q(x) \\)'),
        helpText('\\( \\ \\ \\ \\text{ 2. calculate the ratio } r=\\frac{f(x^*)}{f(x^{(t-1)})}\\frac{q(x^{(t-1)})}{q(x^*)} \\)'),
        helpText('\\( \\ \\ \\ \\text{ 3. set } x^{(t)}=x^* \\text{ with probabilty min(1,}r)\\ ;\\ x^{(t-1)} \\text{ otherwise } \\)'),
        helpText('\\( \\text{Step3 : Check converagence and independence} \\)'))
      })
    
    
    ########  tabPanel 3  ########
    output$plot_cd_para_static = renderPlot({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate({
        react()$df_cd_para %>%
          ggplot(aes(x=time, y=para, group=name, color=name)) +
          geom_line(size=1.5, arrow=arrow()) +
          scale_color_manual(
            values=c('blue', 'red')
          ) +
          ggtitle('parameters converage process') +
          theme_minimal() +
          theme(plot.title=element_text(size=20),
                axis.title=element_text(size=15),
                legend.title=element_text(size=15),
                legend.text=element_text(size=15)) +
          ylab('value')
        })
      })
    
    
    output$plot_cd_loss_static = renderPlot({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      
      isolate({
        react()$df_cd_loss %>%
          ggplot(aes(x=time, y=loss, color=name)) +
          geom_line(size=1.5, arrow=arrow()) +
          scale_color_manual(
            values=c('darkgreen')
          ) +
          ggtitle('loss function') +
          theme_minimal() +
          theme(plot.title=element_text(size=20),
                axis.title=element_text(size=15),
                legend.title=element_text(size=15),
                legend.text=element_text(size=15)) +
          ylab('loss')
        })
      })
    
    sliderValues3_1 = reactive({
      data.frame(
        'mu_hat' = react()$cd_para[length(react()$cd_para[,1]), 1],
        'sigma_hat' = react()$cd_para[length(react()$cd_para[,2]), 2])
      })
    output$values3_1 <- renderTable({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate(sliderValues3_1())
      })
    output$values3_2 <- renderTable({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate(sliderValues3_1())
      })
    
    output$plot_cd_para = renderImage({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate({
        p_anim2 = react()$df_cd_para %>%
          ggplot(aes(x=time, y=para, group=name, color=name)) +
          geom_line(size=1.5) +
          geom_point(size=2) +
          scale_color_manual(
            values=c('blue', 'red')
          ) +
          ggtitle('parameters converage process') +
          theme_minimal() +
          theme(plot.title=element_text(size=20),
                axis.title=element_text(size=15),
                legend.title=element_text(size=15),
                legend.text=element_text(size=15)) +
          ylab('value')
        
        anim2 = p_anim2 +
          transition_reveal(time)
        anim_save("outfile2.gif", animate(anim2, renderer=gifski_renderer(loop=T)))
        # Return a list containing the filename
        list(src="outfile2.gif", contentType="image/gif")
        })
      }, deleteFile = T)
    
    output$plot_cd_loss = renderImage({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate({
        p_anim3 = react()$df_cd_loss %>%
          ggplot(aes(x=time, y=loss, color=name)) +
          geom_line(size=1.5) +
          geom_point(size=2) +
          scale_color_manual(
            values=c('darkgreen')
          ) +
          ggtitle('loss function') +
          theme_minimal() +
          theme(plot.title=element_text(size=20),
                axis.title=element_text(size=15),
                legend.title=element_text(size=15),
                legend.text=element_text(size=15)) +
          ylab('loss')
        
        anim3 = p_anim3 +
          transition_reveal(time)
        anim_save("outfile3.gif", animate(anim3, renderer=gifski_renderer(loop=T)))
        # Return a list containing the filename
        list(src="outfile3.gif", contentType="image/gif")
        })
      }, deleteFile = T)
    
    output$plot_ll_3d = renderPlotly({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate({
        fig1 = plot_ly(x=react()$test_x, y=react()$test_y, z=react()$tmp_m, type='surface') %>%
          layout(scene = list(xaxis=list(title='mu'),
                              yaxis=list(title='sigma'),
                              zaxis=list(title='logL')))
        add_trace(fig1,
                  x=react()$cd_para[length(react()$cd_para[,1]), 1],
                  y=react()$cd_para[length(react()$cd_para[,2]), 2],
                  z=-ll(c(react()$cd_para[length(react()$cd_para[,1]), 1],
                          react()$cd_para[length(react()$cd_para[,2]), 2]),
                        input=react()$acf_x),
                  type='scatter3d', mode='markers', color=I('red')) #add_trace會有Warning, 但可顯示
        })
      })
    
    output$ex3_1 = renderUI({
      withMathJax(helpText('Log-likelihood function : $$\\text{logL}(\\mu, \\sigma)=\\frac{-n}{2}\\text{log}(2\\pi\\sigma^2)+\\sum_{i=1}^{n}\\text{log}(e^{-(x_{i}-\\mu)^2/(2\\sigma^2)}+e^{-(x_{i}+\\mu)^2/(2\\sigma^2)})$$'))
      })
    
    output$plot_ll_heatmap = renderPlotly({
      input$go1_1
      input$go2_1
      input$go2_2
      input$go3_1
      isolate({
        fig2 = plot_ly(x=react()$test_x, y=react()$test_y, z=react()$tmp_m, type = 'heatmap') %>%
          layout(xaxis=list(title='mu'),
                 yaxis=list(title='sigma'))
        add_trace(fig2,
                  x=react()$cd_para[length(react()$cd_para[,1]), 1],
                  y=react()$cd_para[length(react()$cd_para[,2]), 2], 
                  type='scatter', mode='markers', color=I('red')) #add_trace會有Warning, 但可顯示
        })
      })
    
    output$ex3_2 = renderUI({
      withMathJax(
        helpText('$$\\text{logL}(\\mu, \\sigma)=\\frac{-n}{2}\\text{log}(2\\pi\\sigma^2)+\\sum_{i=1}^{n}\\text{log}(e^{-(x_{i}-\\mu)^2/(2\\sigma^2)}+e^{-(x_{i}+\\mu)^2/(2\\sigma^2)})$$'),
        helpText('coodinate descent algorithm : '),
        helpText('\\( \\text{Step1 : Find the negative log-likelihood function, denote } g(\\mu,\\sigma) \\)'),
        helpText('\\( \\text{Step2 : Set initial value of } \\theta^{(0)} = (\\mu^{(0)}, \\sigma^{(0)}) \\)'),
        helpText('\\( \\text{Step3 : Fix } \\sigma^{(t)} \\text{ and get } \\mu^{(t+1)} = arg \\underset{\\mu} min\\ g(\\mu|\\sigma^{(t)}) \\)'),
        helpText('\\( \\text{Step4 : Fix } \\mu^{(t+1)} \\text{ and get } \\sigma^{(t+1)} = arg \\underset{\\sigma} min\\ g(\\sigma|\\mu^{(t+1)}) \\)'),
        helpText('\\( \\text{Step5 : Repeat step3,4 until all of } |\\mu^{(t+1)}-\\mu^{(t)}| < 10^{-10},\\ |\\sigma^{(t+1)}-\\sigma^{(t)}| < 10^{-10} \\)'))
      })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
