:: // Set folder name
set str="ch"

FOR /L %%G IN (1,1,12) DO (
	echo %%G
	IF %%G LSS 10 (			
        mkdir %str%0%%G
    ) ELSE (
        mkdir %str%%%G
    )	 
)
pause


:: // references
:: for /l %x in (1, 1, 100) do (
   :: echo %x
   :: copy %x.txt z:\whatever\etc
:: )