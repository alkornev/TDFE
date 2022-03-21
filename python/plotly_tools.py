def get_sliders(n_frames, dt, fr_duration=100, x_pos=0.0, slider_len=1.0):
    # n_frames= number of frames
    #fr_duration=the duration in milliseconds of each frame
    #x_pos x-coordinate where the slider starts
    #slider_len is a number in (0,1] giving the slider length as a fraction of x-axis length 
    return [dict(steps= [dict(method= 'animate',#Sets the Plotly method to be called when the slider value is changed.
                              args= [ [ 'frame{}'.format(k) ],#Sets the arguments values to be passed to the Plotly,
                                                              #method set in method on slide
                                      dict(mode= 'immediate',
                                           frame= dict( duration=fr_duration, redraw= True ),
                                           transition=dict( duration= 0)
                                          )
                                    ],
                              label='time: {} a.u.'.format(k*dt)
                             ) for k in range(n_frames)], 
                transition= { 'duration': 0 },
                x=x_pos,
                len=slider_len)]

def get_updatemenus(x_pos=0.0, fr_duration=200):
    return [dict(x= x_pos,
                 y= 0,
                 yanchor='top',
                 xanchor= 'right',
                 pad= dict(r= 10, t=40 ),
                 type='buttons',
                 showactive= False,
                 buttons= [dict(label='Play',
                                method='animate',
                                args= [ None,
                                        dict(mode='immediate',
                                             transition= { 'duration': 0 },
                                             fromcurrent= True,
                                             frame= dict( redraw=True, duration=fr_duration)
                                            )
                                       ]
                                ),
                           dict(label='Pause',
                                method='animate',
                                args= [ [None],
                                        dict(mode='immediate',
                                             transition= { 'duration': 0 },
                                             frame= dict( redraw=True, duration=0 )
                                            )
                                       ]
                                )
                           ]
               )
        ]
